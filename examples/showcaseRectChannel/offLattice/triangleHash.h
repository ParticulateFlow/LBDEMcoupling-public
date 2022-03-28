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

#ifndef TRIANGLE_HASH_H
#define TRIANGLE_HASH_H

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiDataField3D.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "particles/particleField3D.h"

namespace plb {

template <typename T>
class TriangleHash {
public:
    TriangleHash(AtomicContainerBlock3D &hashContainer);
    void assignTriangles(TriangularSurfaceMesh<T> const &mesh);
    void bruteReAssignTriangles(TriangularSurfaceMesh<T> const &mesh);
    template <class ParticleFieldT>
    void reAssignTriangles(
        TriangularSurfaceMesh<T> const &mesh, ParticleFieldT &particles,
        std::vector<plint> const &nonParallelVertices);
    void getTriangles(
        Array<T, 2> const &xRange, Array<T, 2> const &yRange, Array<T, 2> const &zRange,
        std::vector<plint> &foundTriangles) const;
    void getTriangles(Box3D const &domain, std::vector<plint> &foundTriangles) const;

private:
    ScalarField3D<std::vector<plint> > &triangles;
    std::vector<Dot3D> &assignedPositions;
};

template <typename T>
class CreateTriangleHash : public BoxProcessingFunctional3D {
public:
    CreateTriangleHash(TriangularSurfaceMesh<T> const &mesh_);
    // Field 0: Hash.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CreateTriangleHash<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    TriangularSurfaceMesh<T> const &mesh;
};

template <typename T, class ParticleFieldT>
class ReAssignTriangleHash : public BoxProcessingFunctional3D {
public:
    // The triangle hash includes the possibilities for triangles
    //   which are too big for proper parallelization, such as, inlet
    //   or outlet triangles. The vertices of these triangles must be
    //   non-parallel, which means, well defined on each processor. The
    //   argument nonParallelVertices must enumerate at least one vertex
    //   of each non-parallel triangle. It is a requirement that all
    //   vertices of all triangles connected to the vertices listed in
    //   nonParallelVertices are non-parallel, i.e. define on each
    //   processor.
    ReAssignTriangleHash(
        TriangularSurfaceMesh<T> const &mesh_, std::vector<plint> const &nonParallelVertices_);
    // Field 0: Hash; Field 1: Particles.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ReAssignTriangleHash<T, ParticleFieldT> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    TriangularSurfaceMesh<T> const &mesh;
    std::vector<plint> nonParallelVertices;
};

template <typename T>
class BruteReAssignTriangleHash : public BoxProcessingFunctional3D {
public:
    BruteReAssignTriangleHash(TriangularSurfaceMesh<T> const &mesh_);
    // Field 0: Hash.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual BruteReAssignTriangleHash<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    TriangularSurfaceMesh<T> const &mesh;
};

}  // namespace plb

#endif  // TRIANGLE_HASH_H
