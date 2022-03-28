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

#ifndef IMMERSED_ADVECTION_DIFFUSION_WALLS_3D_HH
#define IMMERSED_ADVECTION_DIFFUSION_WALLS_3D_HH

#include <vector>

#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "immersedAdvectionDiffusionWalls3D.h"
#include "immersedWalls3D.h"

namespace plb {

/* ******** InamuroAdvectionDiffusionIteration3D ************************************ */

template <typename T, class ScalarFunction>
InamuroAdvectionDiffusionIteration3D<T, ScalarFunction>::InamuroAdvectionDiffusionIteration3D(
    ScalarFunction scalarFunction_, T tau_) :
    scalarFunction(scalarFunction_), tau(tau_)
{ }

template <typename T, class ScalarFunction>
void InamuroAdvectionDiffusionIteration3D<T, ScalarFunction>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ScalarField3D<T> *rhoBar = dynamic_cast<ScalarField3D<T> *>(blocks[0]);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[1]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(container);
    Dot3D location = rhoBar->getLocation();
    ImmersedWallData3D<T> *wallData = dynamic_cast<ImmersedWallData3D<T> *>(container->getData());
    PLB_ASSERT(wallData);

    std::vector<Array<T, 3> > const &vertices = wallData->vertices;
    std::vector<T> const &areas = wallData->areas;
    PLB_ASSERT(vertices.size() == areas.size());
    std::vector<T> deltaS(vertices.size());

    // In this iteration, the source is computed for every vertex.
    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 3> const &vertex = vertices[i];
        Array<plint, 3> intPos(
            (plint)vertex[0] - location.x, (plint)vertex[1] - location.y,
            (plint)vertex[2] - location.z);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> zLim(
            (vertex[2] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 3> fracPos(
            util::frac(vertex[0]), util::frac(vertex[1]), util::frac(vertex[2]));
        T averageRhoBar = T();
        // Use the weighting function to compute the average density
        // on the surface vertex.
        // x   x . x   x
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                    Array<T, 3> r((T)dx - fracPos[0], (T)dy - fracPos[1], (T)dz - fracPos[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageRhoBar += W * nextRhoBar;
                }
            }
        }
        T wallScalar = scalarFunction(vertex);
        deltaS[i] = areas[i] * ((wallScalar - (T)1.0) - averageRhoBar);
    }

    // In this iteration, the source is applied from every vertex to the grid nodes.
    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 3> const &vertex = vertices[i];
        Array<plint, 3> intPos(
            (plint)vertex[0] - location.x, (plint)vertex[1] - location.y,
            (plint)vertex[2] - location.z);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> zLim(
            (vertex[2] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 3> fracPos(
            util::frac(vertex[0]), util::frac(vertex[1]), util::frac(vertex[2]));
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                    Array<T, 3> r((T)dx - fracPos[0], (T)dy - fracPos[1], (T)dz - fracPos[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    nextRhoBar += tau * W * deltaS[i];
                    rhoBar->get(pos[0], pos[1], pos[2]) = nextRhoBar;
                }
            }
        }
    }
}

template <typename T, class ScalarFunction>
InamuroAdvectionDiffusionIteration3D<T, ScalarFunction>
    *InamuroAdvectionDiffusionIteration3D<T, ScalarFunction>::clone() const
{
    return new InamuroAdvectionDiffusionIteration3D<T, ScalarFunction>(*this);
}

template <typename T, class ScalarFunction>
void InamuroAdvectionDiffusionIteration3D<T, ScalarFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // RhoBar
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class ScalarFunction>
BlockDomain::DomainT InamuroAdvectionDiffusionIteration3D<T, ScalarFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** IndexedInamuroAdvectionDiffusionIteration3D ************************************ */

template <typename T, class ScalarFunction>
IndexedInamuroAdvectionDiffusionIteration3D<T, ScalarFunction>::
    IndexedInamuroAdvectionDiffusionIteration3D(ScalarFunction scalarFunction_, T tau_) :
    scalarFunction(scalarFunction_), tau(tau_)
{ }

template <typename T, class ScalarFunction>
void IndexedInamuroAdvectionDiffusionIteration3D<T, ScalarFunction>::processGenericBlocks(
    [[maybe_unused]] Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ScalarField3D<T> *rhoBar = dynamic_cast<ScalarField3D<T> *>(blocks[0]);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[1]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(container);
    Dot3D location = rhoBar->getLocation();
    ImmersedWallData3D<T> *wallData = dynamic_cast<ImmersedWallData3D<T> *>(container->getData());
    PLB_ASSERT(wallData);

    std::vector<Array<T, 3> > const &vertices = wallData->vertices;
    std::vector<T> const &areas = wallData->areas;
    PLB_ASSERT(vertices.size() == areas.size());
    std::vector<T> deltaS(vertices.size());
    std::vector<pluint> const &globalVertexIds = wallData->globalVertexIds;
    PLB_ASSERT(vertices.size() == globalVertexIds.size());

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 3> const &vertex = vertices[i];
        Array<plint, 3> intPos(
            (plint)vertex[0] - location.x, (plint)vertex[1] - location.y,
            (plint)vertex[2] - location.z);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> zLim(
            (vertex[2] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 3> fracPos(
            util::frac(vertex[0]), util::frac(vertex[1]), util::frac(vertex[2]));
        T averageRhoBar = T();
        // x   x . x   x
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                    Array<T, 3> r((T)dx - fracPos[0], (T)dy - fracPos[1], (T)dz - fracPos[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageRhoBar += W * nextRhoBar;
                }
            }
        }
        T wallScalar = scalarFunction(globalVertexIds[i]);
        deltaS[i] = areas[i] * ((wallScalar - (T)1.0) - averageRhoBar);
    }

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 3> const &vertex = vertices[i];
        Array<plint, 3> intPos(
            (plint)vertex[0] - location.x, (plint)vertex[1] - location.y,
            (plint)vertex[2] - location.z);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> zLim(
            (vertex[2] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 3> fracPos(
            util::frac(vertex[0]), util::frac(vertex[1]), util::frac(vertex[2]));
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                    Array<T, 3> r((T)dx - fracPos[0], (T)dy - fracPos[1], (T)dz - fracPos[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    nextRhoBar += tau * W * deltaS[i];
                    rhoBar->get(pos[0], pos[1], pos[2]) = nextRhoBar;
                }
            }
        }
    }
}

template <typename T, class ScalarFunction>
IndexedInamuroAdvectionDiffusionIteration3D<T, ScalarFunction>
    *IndexedInamuroAdvectionDiffusionIteration3D<T, ScalarFunction>::clone() const
{
    return new IndexedInamuroAdvectionDiffusionIteration3D<T, ScalarFunction>(*this);
}

template <typename T, class ScalarFunction>
void IndexedInamuroAdvectionDiffusionIteration3D<T, ScalarFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // RhoBar
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class ScalarFunction>
BlockDomain::DomainT IndexedInamuroAdvectionDiffusionIteration3D<T, ScalarFunction>::appliesTo()
    const
{
    return BlockDomain::bulk;
}

/* ******** ConstScalarInamuroAdvectionDiffusionIteration3D ************************************ */

template <typename T>
ConstScalarInamuroAdvectionDiffusionIteration3D<T>::ConstScalarInamuroAdvectionDiffusionIteration3D(
    T wallScalar_, T tau_) :
    wallScalar(wallScalar_), tau(tau_)
{ }

template <typename T>
void ConstScalarInamuroAdvectionDiffusionIteration3D<T>::processGenericBlocks(
    [[maybe_unused]] Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ScalarField3D<T> *rhoBar = dynamic_cast<ScalarField3D<T> *>(blocks[0]);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[1]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(container);
    Dot3D location = rhoBar->getLocation();
    ImmersedWallData3D<T> *wallData = dynamic_cast<ImmersedWallData3D<T> *>(container->getData());
    PLB_ASSERT(wallData);
    std::vector<Array<T, 3> > const &vertices = wallData->vertices;
    std::vector<T> const &areas = wallData->areas;
    PLB_ASSERT(vertices.size() == areas.size());
    std::vector<T> deltaS(vertices.size());

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 3> const &vertex = vertices[i];
        Array<plint, 3> intPos(
            (plint)vertex[0] - location.x, (plint)vertex[1] - location.y,
            (plint)vertex[2] - location.z);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> zLim(
            (vertex[2] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 3> fracPos(
            util::frac(vertex[0]), util::frac(vertex[1]), util::frac(vertex[2]));
        T averageRhoBar = T();
        // x   x . x   x
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                    Array<T, 3> r((T)dx - fracPos[0], (T)dy - fracPos[1], (T)dz - fracPos[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageRhoBar += W * nextRhoBar;
                }
            }
        }
        deltaS[i] = areas[i] * ((wallScalar - (T)1.0) - averageRhoBar);
    }

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 3> const &vertex = vertices[i];
        Array<plint, 3> intPos(
            (plint)vertex[0] - location.x, (plint)vertex[1] - location.y,
            (plint)vertex[2] - location.z);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> zLim(
            (vertex[2] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 3> fracPos(
            util::frac(vertex[0]), util::frac(vertex[1]), util::frac(vertex[2]));
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                    Array<T, 3> r((T)dx - fracPos[0], (T)dy - fracPos[1], (T)dz - fracPos[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    nextRhoBar += tau * W * deltaS[i];
                    rhoBar->get(pos[0], pos[1], pos[2]) = nextRhoBar;
                }
            }
        }
    }
}

template <typename T>
ConstScalarInamuroAdvectionDiffusionIteration3D<T>
    *ConstScalarInamuroAdvectionDiffusionIteration3D<T>::clone() const
{
    return new ConstScalarInamuroAdvectionDiffusionIteration3D<T>(*this);
}

template <typename T>
void ConstScalarInamuroAdvectionDiffusionIteration3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // RhoBar
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT ConstScalarInamuroAdvectionDiffusionIteration3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** InstantiateImmersedAdvectionDiffusionWallData3D ************************************ */

template <typename T>
InstantiateImmersedAdvectionDiffusionWallData3D<T>::InstantiateImmersedAdvectionDiffusionWallData3D(
    std::vector<Array<T, 3> > const &vertices_, std::vector<T> const &areas_,
    std::vector<Array<T, 3> > const &normals_) :
    vertices(vertices_), areas(areas_), normals(normals_)
{
    PLB_ASSERT(vertices.size() == areas.size());
    PLB_ASSERT(normals.size() == 0 || normals.size() == areas.size());
}

template <typename T>
void InstantiateImmersedAdvectionDiffusionWallData3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[0]);
    PLB_ASSERT(container);
    bool useNormals = normals.size() > 0;
    Dot3D location = container->getLocation();
    ImmersedWallData3D<T> *wallData = new ImmersedWallData3D<T>;
    Box3D extendedEnvelope(domain.enlarge(2).shift(location.x, location.y, location.z));

    for (pluint i = 0; i < vertices.size(); ++i) {
        if (contained(vertices[i], extendedEnvelope)) {
            wallData->vertices.push_back(vertices[i]);
            wallData->areas.push_back(areas[i]);
            if (useNormals) {
                wallData->normals.push_back(normals[i]);
            }
            wallData->globalVertexIds.push_back(i);
        }
    }
    wallData->g.clear();
    wallData->flags = std::vector<int>(wallData->vertices.size(), 0);
    container->setData(wallData);
}

template <typename T>
InstantiateImmersedAdvectionDiffusionWallData3D<T>
    *InstantiateImmersedAdvectionDiffusionWallData3D<T>::clone() const
{
    return new InstantiateImmersedAdvectionDiffusionWallData3D<T>(*this);
}

template <typename T>
void InstantiateImmersedAdvectionDiffusionWallData3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT InstantiateImmersedAdvectionDiffusionWallData3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D
 * ************************************ */

template <typename T>
InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D<T>::
    InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D(
        std::vector<Array<T, 3> > const &vertices_, std::vector<T> const &areas_,
        int solutionDomainFlag_) :
    vertices(vertices_), areas(areas_), solutionDomainFlag(solutionDomainFlag_)
{
    PLB_ASSERT(vertices.size() == areas.size());
}

template <typename T>
void InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[0]);
    PLB_ASSERT(container);
    Dot3D location = container->getLocation();

    ScalarField3D<int> *flagMatrix = dynamic_cast<ScalarField3D<int> *>(blocks[1]);
    PLB_ASSERT(flagMatrix);
    Dot3D ofsFlag = computeRelativeDisplacement(*container, *flagMatrix);
    Array<plint, 3> flagDispl(ofsFlag.x, ofsFlag.y, ofsFlag.z);

    ImmersedWallData3D<T> *wallData = new ImmersedWallData3D<T>;
    Box3D extendedEnvelope(domain.enlarge(2).shift(location.x, location.y, location.z));

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 3> vertex = vertices[i];
        if (contained(vertex, extendedEnvelope)) {
            wallData->vertices.push_back(vertex);
            wallData->areas.push_back(areas[i]);
            wallData->globalVertexIds.push_back(i);
            Array<plint, 3> intPos(
                (plint)vertex[0] - location.x, (plint)vertex[1] - location.y,
                (plint)vertex[2] - location.z);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> zLim(
                (vertex[2] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            bool hasSolutionDomainNeighbor = false;
            for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
                for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                    for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                        Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz) + flagDispl);
                        if (flagMatrix->get(pos[0], pos[1], pos[2]) == solutionDomainFlag) {
                            hasSolutionDomainNeighbor = true;
                        }
                    }
                }
            }
            if (hasSolutionDomainNeighbor) {
                wallData->flags.push_back(0);
            } else {
                wallData->flags.push_back(1);
            }
        }
    }
    wallData->g.clear();
    container->setData(wallData);
}

template <typename T>
InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D<T>
    *InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D<T>::clone() const
{
    return new InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D<T>(*this);
}

template <typename T>
void InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
    modified[1] = modif::nothing;          // Flag matrix.
}

template <typename T>
BlockDomain::DomainT InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D<T>::appliesTo()
    const
{
    return BlockDomain::bulk;
}

/* ******** InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D
 * ************************************ */

template <typename T>
InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D<T>::
    InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D(
        std::vector<Array<T, 3> > const &vertices_, std::vector<T> const &areas_,
        std::vector<int> const &flags_) :
    vertices(vertices_), areas(areas_), flags(flags_)
{
    PLB_ASSERT(vertices.size() == areas.size());
    PLB_ASSERT(vertices.size() == flags.size());
}

template <typename T>
void InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[0]);
    PLB_ASSERT(container);
    Dot3D location = container->getLocation();
    ImmersedWallData3D<T> *wallData = new ImmersedWallData3D<T>;
    Box3D extendedEnvelope(domain.enlarge(2).shift(location.x, location.y, location.z));

    for (pluint i = 0; i < vertices.size(); ++i) {
        if (contained(vertices[i], extendedEnvelope)) {
            wallData->vertices.push_back(vertices[i]);
            wallData->areas.push_back(areas[i]);
            wallData->flags.push_back(flags[i]);
            wallData->globalVertexIds.push_back(i);
        }
    }
    wallData->g.clear();
    container->setData(wallData);
}

template <typename T>
InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D<T>
    *InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D<T>::clone() const
{
    return new InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D<T>(*this);
}

template <typename T>
void InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT
    InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** BoxAdvectionDiffusionRhoBarJfunctional3D ************************************ */

template <
    typename T, template <typename U> class Descriptor, template <typename V> class AD_Descriptor>
void BoxAdvectionDiffusionRhoBarJfunctional3D<T, Descriptor, AD_Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_ASSERT(fields.size() == 4);

    BlockLattice3D<T, AD_Descriptor> *ad_lattice =
        dynamic_cast<BlockLattice3D<T, AD_Descriptor> *>(fields[0]);
    PLB_ASSERT(ad_lattice);
    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(fields[1]);
    PLB_ASSERT(lattice);
    ScalarField3D<T> *rhoBar = dynamic_cast<ScalarField3D<T> *>(fields[2]);
    PLB_ASSERT(rhoBar);
    TensorField3D<T, 3> *j = dynamic_cast<TensorField3D<T, 3> *>(fields[3]);
    PLB_ASSERT(j);

    Dot3D ofsL = computeRelativeDisplacement(*ad_lattice, *lattice);
    Dot3D ofsR = computeRelativeDisplacement(*ad_lattice, *rhoBar);
    Dot3D ofsJ = computeRelativeDisplacement(*ad_lattice, *j);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T rho = ad_lattice->get(iX, iY, iZ).computeDensity();
                Array<T, 3> u;
                lattice->get(iX + ofsL.x, iY + ofsL.y, iZ + ofsL.z).computeVelocity(u);

                rhoBar->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z) = AD_Descriptor<T>::rhoBar(rho);
                j->get(iX + ofsJ.x, iY + ofsJ.y, iZ + ofsJ.z) = rho * u;
            }
        }
    }
}

template <
    typename T, template <typename U> class Descriptor, template <typename V> class AD_Descriptor>
BoxAdvectionDiffusionRhoBarJfunctional3D<T, Descriptor, AD_Descriptor>
    *BoxAdvectionDiffusionRhoBarJfunctional3D<T, Descriptor, AD_Descriptor>::clone() const
{
    return new BoxAdvectionDiffusionRhoBarJfunctional3D<T, Descriptor, AD_Descriptor>(*this);
}

template <
    typename T, template <typename U> class Descriptor, template <typename V> class AD_Descriptor>
void BoxAdvectionDiffusionRhoBarJfunctional3D<T, Descriptor, AD_Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Advection-diffusion lattice.
    modified[1] = modif::nothing;          // Fluid lattice.
    modified[2] = modif::staticVariables;  // rhoBar
    modified[3] = modif::staticVariables;  // j
}

}  // namespace plb

#endif  // IMMERSED_ADVECTION_DIFFUSION_WALLS_3D_HH
