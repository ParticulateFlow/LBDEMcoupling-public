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

#ifndef TRIANGLE_HASH_HH
#define TRIANGLE_HASH_HH

#include <algorithm>

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "core/globalDefs.h"
#include "offLattice/triangleHash.h"

namespace plb {

struct TriangleHashData : public ContainerBlockData {
    TriangleHashData(plint nx, plint ny, plint nz, Dot3D const &location) : triangles(nx, ny, nz)
    {
        triangles.setLocation(location);
    }
    virtual TriangleHashData *clone() const
    {
        return new TriangleHashData(*this);
    }
    ScalarField3D<std::vector<plint> > triangles;
    std::vector<Dot3D> assignedPositions;
};

/* ******** class TriangleHash ********************************************* */

template <typename T>
TriangleHash<T>::TriangleHash(AtomicContainerBlock3D &hashContainer) :
    triangles(dynamic_cast<TriangleHashData *>(hashContainer.getData())->triangles),
    assignedPositions(dynamic_cast<TriangleHashData *>(hashContainer.getData())->assignedPositions)
{ }

template <typename T>
void TriangleHash<T>::getTriangles(
    Array<T, 2> const &xRange, Array<T, 2> const &yRange, Array<T, 2> const &zRange,
    std::vector<plint> &foundTriangles) const
{
    // Fit onto the grid by making it bigger, to be sure the triangle
    //   is never missed through round-off errors.
    Box3D discreteRange(
        (plint)xRange[0], (plint)xRange[1] + 1, (plint)yRange[0], (plint)yRange[1] + 1,
        (plint)zRange[0], (plint)zRange[1] + 1);
    getTriangles(discreteRange, foundTriangles);
}

template <typename T>
void TriangleHash<T>::getTriangles(Box3D const &domain, std::vector<plint> &foundTriangles) const
{
    Dot3D location(triangles.getLocation());
    // Convert to local coordinates.
    Box3D shifted(domain.shift(-location.x, -location.y, -location.z));
    foundTriangles.clear();
    Box3D inters;
    if (intersect(shifted, triangles.getBoundingBox(), inters)) {
        for (plint iX = inters.x0; iX <= inters.x1; ++iX) {
            for (plint iY = inters.y0; iY <= inters.y1; ++iY) {
                for (plint iZ = inters.z0; iZ <= inters.z1; ++iZ) {
                    std::vector<plint> const &newTriangles = triangles.get(iX, iY, iZ);
                    foundTriangles.insert(
                        foundTriangles.end(), newTriangles.begin(), newTriangles.end());
                }
            }
        }
        std::sort(foundTriangles.begin(), foundTriangles.end());
        foundTriangles.erase(
            unique(foundTriangles.begin(), foundTriangles.end()), foundTriangles.end());
    }
}

template <typename T>
void TriangleHash<T>::assignTriangles(TriangularSurfaceMesh<T> const &mesh)
{
    Dot3D location(triangles.getLocation());
    assignedPositions.clear();
    for (plint iTriangle = 0; iTriangle < mesh.getNumTriangles(); ++iTriangle) {
        Array<T, 3> const &vertex0 = mesh.getVertex(iTriangle, 0);
        Array<T, 3> const &vertex1 = mesh.getVertex(iTriangle, 1);
        Array<T, 3> const &vertex2 = mesh.getVertex(iTriangle, 2);

        Array<T, 2> xRange(
            std::min(vertex0[0], std::min(vertex1[0], vertex2[0])),
            std::max(vertex0[0], std::max(vertex1[0], vertex2[0])));
        Array<T, 2> yRange(
            std::min(vertex0[1], std::min(vertex1[1], vertex2[1])),
            std::max(vertex0[1], std::max(vertex1[1], vertex2[1])));
        Array<T, 2> zRange(
            std::min(vertex0[2], std::min(vertex1[2], vertex2[2])),
            std::max(vertex0[2], std::max(vertex1[2], vertex2[2])));

        // Fit onto the grid by making it bigger, to be sure the triangle
        //   is never missed through round-off errors.
        Box3D discreteRange(
            (plint)xRange[0], (plint)xRange[1] + 1, (plint)yRange[0], (plint)yRange[1] + 1,
            (plint)zRange[0], (plint)zRange[1] + 1);
        // Convert to local coordinates.
        discreteRange = discreteRange.shift(-location.x, -location.y, -location.z);
        Box3D inters;
        if (intersect(discreteRange, triangles.getBoundingBox(), inters)) {
            for (plint iX = inters.x0; iX <= inters.x1; ++iX) {
                for (plint iY = inters.y0; iY <= inters.y1; ++iY) {
                    for (plint iZ = inters.z0; iZ <= inters.z1; ++iZ) {
                        if (triangles.get(iX, iY, iZ).empty()) {
                            assignedPositions.push_back(Dot3D(iX, iY, iZ));
                        }
                        triangles.get(iX, iY, iZ).push_back(iTriangle);
                    }
                }
            }
        }
    }
}

template <typename T>
template <class ParticleFieldT>
void TriangleHash<T>::reAssignTriangles(
    TriangularSurfaceMesh<T> const &mesh, ParticleFieldT &particles,
    std::vector<plint> const &nonParallelVertices)
{
    // Create domain from which particles are going to be retrieved.
    Box3D domain(triangles.getBoundingBox());
    // First of all, remove old triangles.
    /*
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                triangles.get(iX,iY,iZ).clear();
            }
        }
    }
    */
    for (pluint iAssigned = 0; iAssigned < assignedPositions.size(); ++iAssigned) {
        Dot3D pos(assignedPositions[iAssigned]);
        triangles.get(pos.x, pos.y, pos.z).clear();
    }
    assignedPositions.clear();
    // Particles have a bigger envelope than triangles (which have
    //   envelope with 2 for Guo for example). The domain must now
    //   be translated into the local coordinates of the particles.
    Dot3D offset = computeRelativeDisplacement(triangles, particles);
    domain = domain.shift(offset.x, offset.y, offset.z);
    // Enlarge by one cell, because triangles belong to the hash of
    //   a given AtomicBlock even when one of their vertices is out-
    //   side the AtomicBlock by maximally one cell.
    domain.enlarge(1);
    std::vector<typename ParticleFieldT::ParticleT *> found;
    particles.findParticles(domain, found);
    std::set<plint> triangleIds;
    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        plint vertexId = found[iParticle]->getTag();
        std::vector<plint> newTriangles(mesh.getNeighborTriangleIds(vertexId));
        triangleIds.insert(newTriangles.begin(), newTriangles.end());
    }
    for (pluint iVertex = 0; iVertex < nonParallelVertices.size(); ++iVertex) {
        plint vertexId = nonParallelVertices[iVertex];
        std::vector<plint> newTriangles(mesh.getNeighborTriangleIds(vertexId));
        triangleIds.insert(newTriangles.begin(), newTriangles.end());
    }

    Dot3D location(triangles.getLocation());
    std::set<plint>::const_iterator it = triangleIds.begin();
    for (; it != triangleIds.end(); ++it) {
        plint iTriangle = *it;
        Array<T, 3> const &vertex0 = mesh.getVertex(iTriangle, 0);
        Array<T, 3> const &vertex1 = mesh.getVertex(iTriangle, 1);
        Array<T, 3> const &vertex2 = mesh.getVertex(iTriangle, 2);

        Array<T, 2> xRange(
            std::min(vertex0[0], std::min(vertex1[0], vertex2[0])),
            std::max(vertex0[0], std::max(vertex1[0], vertex2[0])));
        Array<T, 2> yRange(
            std::min(vertex0[1], std::min(vertex1[1], vertex2[1])),
            std::max(vertex0[1], std::max(vertex1[1], vertex2[1])));
        Array<T, 2> zRange(
            std::min(vertex0[2], std::min(vertex1[2], vertex2[2])),
            std::max(vertex0[2], std::max(vertex1[2], vertex2[2])));

        // Fit onto the grid by making it bigger, to be sure the triangle
        //   is never missed through round-off errors.
        Box3D discreteRange(
            (plint)xRange[0], (plint)xRange[1] + 1, (plint)yRange[0], (plint)yRange[1] + 1,
            (plint)zRange[0], (plint)zRange[1] + 1);
        // Convert to local coordinates.
        discreteRange = discreteRange.shift(-location.x, -location.y, -location.z);
        Box3D inters;
        if (intersect(discreteRange, triangles.getBoundingBox(), inters)) {
            for (plint iX = inters.x0; iX <= inters.x1; ++iX) {
                for (plint iY = inters.y0; iY <= inters.y1; ++iY) {
                    for (plint iZ = inters.z0; iZ <= inters.z1; ++iZ) {
                        if (triangles.get(iX, iY, iZ).empty()) {
                            assignedPositions.push_back(Dot3D(iX, iY, iZ));
                        }
                        triangles.get(iX, iY, iZ).push_back(iTriangle);
                    }
                }
            }
        }
    }
}

template <typename T>
void TriangleHash<T>::bruteReAssignTriangles(TriangularSurfaceMesh<T> const &mesh)
{
    PLB_ASSERT(false);
    // First of all, remove old triangles.
    Box3D domain(triangles.getBoundingBox());
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                triangles.get(iX, iY, iZ).clear();
            }
        }
    }
    Dot3D location(triangles.getLocation());
    for (plint iTriangle = 0; iTriangle < mesh.getNumTriangles(); ++iTriangle) {
        if (mesh.isValidVertex(iTriangle, 0) && mesh.isValidVertex(iTriangle, 1)
            && mesh.isValidVertex(iTriangle, 2))
        {
            Array<T, 3> const &vertex0 = mesh.getVertex(iTriangle, 0);
            Array<T, 3> const &vertex1 = mesh.getVertex(iTriangle, 1);
            Array<T, 3> const &vertex2 = mesh.getVertex(iTriangle, 2);

            Array<T, 2> xRange(
                std::min(vertex0[0], std::min(vertex1[0], vertex2[0])),
                std::max(vertex0[0], std::max(vertex1[0], vertex2[0])));
            Array<T, 2> yRange(
                std::min(vertex0[1], std::min(vertex1[1], vertex2[1])),
                std::max(vertex0[1], std::max(vertex1[1], vertex2[1])));
            Array<T, 2> zRange(
                std::min(vertex0[2], std::min(vertex1[2], vertex2[2])),
                std::max(vertex0[2], std::max(vertex1[2], vertex2[2])));

            // Fit onto the grid by making it bigger, to be sure the triangle
            //   is never missed through round-off errors.
            Box3D discreteRange(
                (plint)xRange[0], (plint)xRange[1] + 1, (plint)yRange[0], (plint)yRange[1] + 1,
                (plint)zRange[0], (plint)zRange[1] + 1);
            // Convert to local coordinates.
            discreteRange = discreteRange.shift(-location.x, -location.y, -location.z);
            Box3D inters;
            if (intersect(discreteRange, triangles.getBoundingBox(), inters)) {
                for (plint iX = inters.x0; iX <= inters.x1; ++iX) {
                    for (plint iY = inters.y0; iY <= inters.y1; ++iY) {
                        for (plint iZ = inters.z0; iZ <= inters.z1; ++iZ) {
                            triangles.get(iX, iY, iZ).push_back(iTriangle);
                        }
                    }
                }
            }
        }
    }
}

/* ******** CreateTriangleHash ************************************ */

template <typename T>
CreateTriangleHash<T>::CreateTriangleHash(TriangularSurfaceMesh<T> const &mesh_) : mesh(mesh_)
{ }

template <typename T>
void CreateTriangleHash<T>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[0]);
    PLB_ASSERT(container);
    TriangleHashData *hashData = new TriangleHashData(
        container->getNx(), container->getNy(), container->getNz(), container->getLocation());
    container->setData(hashData);
    TriangleHash<T>(*container).assignTriangles(mesh);
}

template <typename T>
CreateTriangleHash<T> *CreateTriangleHash<T>::clone() const
{
    return new CreateTriangleHash<T>(*this);
}

template <typename T>
void CreateTriangleHash<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with hash data.
}

template <typename T>
BlockDomain::DomainT CreateTriangleHash<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ReAssignTriangleHash ************************************ */

template <typename T, class ParticleFieldT>
ReAssignTriangleHash<T, ParticleFieldT>::ReAssignTriangleHash(
    TriangularSurfaceMesh<T> const &mesh_, std::vector<plint> const &nonParallelVertices_) :
    mesh(mesh_), nonParallelVertices(nonParallelVertices_)
{ }

template <typename T, class ParticleFieldT>
void ReAssignTriangleHash<T, ParticleFieldT>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[0]);
    PLB_ASSERT(container);
    ParticleFieldT *particles = dynamic_cast<ParticleFieldT *>(blocks[1]);
    PLB_ASSERT(particles);
    TriangleHash<T>(*container).reAssignTriangles(mesh, *particles, nonParallelVertices);
}

template <typename T, class ParticleFieldT>
ReAssignTriangleHash<T, ParticleFieldT> *ReAssignTriangleHash<T, ParticleFieldT>::clone() const
{
    return new ReAssignTriangleHash<T, ParticleFieldT>(*this);
}

template <typename T, class ParticleFieldT>
void ReAssignTriangleHash<T, ParticleFieldT>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with hash data.
    modified[1] = modif::nothing;          // Vertex-Particles.
}

template <typename T, class ParticleFieldT>
BlockDomain::DomainT ReAssignTriangleHash<T, ParticleFieldT>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** BruteReAssignTriangleHash ************************************ */

template <typename T>
BruteReAssignTriangleHash<T>::BruteReAssignTriangleHash(TriangularSurfaceMesh<T> const &mesh_) :
    mesh(mesh_)
{ }

template <typename T>
void BruteReAssignTriangleHash<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(blocks[0]);
    PLB_ASSERT(container);
    TriangleHash<T>(*container).bruteReAssignTriangles(mesh);
}

template <typename T>
BruteReAssignTriangleHash<T> *BruteReAssignTriangleHash<T>::clone() const
{
    return new BruteReAssignTriangleHash<T>(*this);
}

template <typename T>
void BruteReAssignTriangleHash<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with hash data.
    modified[1] = modif::nothing;          // Vertex-Particles.
}

template <typename T>
BlockDomain::DomainT BruteReAssignTriangleHash<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // TRIANGLE_HASH_HH
