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

#ifndef VISUAL_PARTICLE_FUNCTIONAL_3D_HH
#define VISUAL_PARTICLE_FUNCTIONAL_3D_HH

#include <algorithm>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "finiteDifference/interpolations2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/immersedWalls3D.h"
#include "offLattice/voxelizer.h"
#include "particles/particleField3D.h"
#include "particles/visualParticleFunctional3D.h"
#include "sitmo/prng_engine.hpp"

namespace plb {

/* ******** CreateParticleFromVertex3D *********************************** */

template <typename T, template <typename U> class Descriptor, class ParticleT>
CreateParticleFromVertex3D<T, Descriptor, ParticleT>::CreateParticleFromVertex3D(
    TriangularSurfaceMesh<T> const &mesh_) :
    mesh(mesh_)
{ }

template <typename T, template <typename U> class Descriptor, class ParticleT>
void CreateParticleFromVertex3D<T, Descriptor, ParticleT>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField3D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);

    for (plint iVertex = 0; iVertex < mesh.getNumVertices(); ++iVertex) {
        Array<T, 3> vertex(mesh.getVertex(iVertex));
        ParticleT *particle = new ParticleT(iVertex, vertex);
        particleField.addParticle(domain, particle);
    }
}

template <typename T, template <typename U> class Descriptor, class ParticleT>
CreateParticleFromVertex3D<T, Descriptor, ParticleT>
    *CreateParticleFromVertex3D<T, Descriptor, ParticleT>::clone() const
{
    return new CreateParticleFromVertex3D<T, Descriptor, ParticleT>(*this);
}

template <typename T, template <typename U> class Descriptor, class ParticleT>
BlockDomain::DomainT CreateParticleFromVertex3D<T, Descriptor, ParticleT>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, class ParticleT>
void CreateParticleFromVertex3D<T, Descriptor, ParticleT>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** CreateParticleFromVertexNoTagging3D *********************************** */

template <typename T, template <typename U> class Descriptor>
CreateParticleFromVertexNoTagging3D<T, Descriptor>::CreateParticleFromVertexNoTagging3D(
    TriangularSurfaceMesh<T> const &mesh_, Particle3D<T, Descriptor> *particleTemplate_) :
    mesh(&mesh_), particleTemplate(particleTemplate_)
{ }

template <typename T, template <typename U> class Descriptor>
CreateParticleFromVertexNoTagging3D<T, Descriptor>::~CreateParticleFromVertexNoTagging3D()
{
    delete particleTemplate;
}

template <typename T, template <typename U> class Descriptor>
CreateParticleFromVertexNoTagging3D<T, Descriptor>::CreateParticleFromVertexNoTagging3D(
    CreateParticleFromVertexNoTagging3D<T, Descriptor> const &rhs) :
    mesh(rhs.mesh), particleTemplate(rhs.particleTemplate->clone())
{ }

template <typename T, template <typename U> class Descriptor>
CreateParticleFromVertexNoTagging3D<T, Descriptor>
    &CreateParticleFromVertexNoTagging3D<T, Descriptor>::operator=(
        CreateParticleFromVertexNoTagging3D<T, Descriptor> const &rhs)
{
    CreateParticleFromVertexNoTagging3D<T, Descriptor>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void CreateParticleFromVertexNoTagging3D<T, Descriptor>::swap(
    CreateParticleFromVertexNoTagging3D<T, Descriptor> &rhs)
{
    std::swap(mesh, rhs.mesh);
    std::swap(particleTemplate, rhs.particleTemplate);
}

template <typename T, template <typename U> class Descriptor>
void CreateParticleFromVertexNoTagging3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField3D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);

    for (plint iVertex = 0; iVertex < mesh->getNumVertices(); ++iVertex) {
        Array<T, 3> vertex(mesh->getVertex(iVertex));
        Particle3D<T, Descriptor> *newParticle = particleTemplate->clone();
        newParticle->getPosition() = vertex;
        particleField.addParticle(domain, newParticle);
    }
}

template <typename T, template <typename U> class Descriptor>
CreateParticleFromVertexNoTagging3D<T, Descriptor>
    *CreateParticleFromVertexNoTagging3D<T, Descriptor>::clone() const
{
    return new CreateParticleFromVertexNoTagging3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CreateParticleFromVertexNoTagging3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** CopyVertexToParticle3D *********************************** */

template <typename T, template <typename U> class Descriptor, class ParticleT>
CopyVertexToParticle3D<T, Descriptor, ParticleT>::CopyVertexToParticle3D(
    TriangularSurfaceMesh<T> const &mesh_) :
    mesh(mesh_)
{ }

template <typename T, template <typename U> class Descriptor, class ParticleT>
void CopyVertexToParticle3D<T, Descriptor, ParticleT>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField3D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);

    std::vector<Particle3D<T, Descriptor> *> found;
    particleField.findParticles(domain, found);
    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        Particle3D<T, Descriptor> &particle = found[iParticle];
        plint tag = particle.getTag();
        Array<T, 3> const &vertexPos = mesh.getVertex(tag);
        particle.getPosition() = vertexPos;
    }
}

template <typename T, template <typename U> class Descriptor, class ParticleT>
CopyVertexToParticle3D<T, Descriptor, ParticleT>
    *CopyVertexToParticle3D<T, Descriptor, ParticleT>::clone() const
{
    return new CopyVertexToParticle3D<T, Descriptor, ParticleT>(*this);
}

template <typename T, template <typename U> class Descriptor, class ParticleT>
BlockDomain::DomainT CopyVertexToParticle3D<T, Descriptor, ParticleT>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor, class ParticleT>
void CopyVertexToParticle3D<T, Descriptor, ParticleT>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** CloneParticles3D *********************************** */

template <typename T, template <typename U> class Descriptor>
CloneParticles3D<T, Descriptor>::CloneParticles3D(Particle3D<T, Descriptor> *particleTemplate_) :
    particleTemplate(particleTemplate_)
{ }

template <typename T, template <typename U> class Descriptor>
CloneParticles3D<T, Descriptor>::CloneParticles3D(CloneParticles3D<T, Descriptor> const &rhs) :
    particleTemplate(rhs.particleTemplate->clone())
{ }

template <typename T, template <typename U> class Descriptor>
CloneParticles3D<T, Descriptor>::~CloneParticles3D()
{
    delete particleTemplate;
}

template <typename T, template <typename U> class Descriptor>
CloneParticles3D<T, Descriptor> &CloneParticles3D<T, Descriptor>::operator=(
    CloneParticles3D<T, Descriptor> const &rhs)
{
    CloneParticles3D<T, Descriptor>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void CloneParticles3D<T, Descriptor>::swap(CloneParticles3D<T, Descriptor> &rhs)
{
    std::swap(particleTemplate, rhs.particleTemplate);
}

template <typename T, template <typename U> class Descriptor>
void CloneParticles3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField3D<T, Descriptor> &originalField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    ParticleField3D<T, Descriptor> &clonedField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[1]);

    std::vector<Particle3D<T, Descriptor> *> found;
    originalField.findParticles(domain, found);

    Dot3D offset = computeRelativeDisplacement(originalField, clonedField);
    Box3D clonedDomain(domain.shift(offset.x, offset.y, offset.z));

    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        Particle3D<T, Descriptor> &originalParticle = *found[iParticle];
        Particle3D<T, Descriptor> *clonedParticle = particleTemplate->clone();
        clonedParticle->reset(originalParticle.getPosition());
        clonedParticle->setTag(originalParticle.getTag());
        clonedField.addParticle(clonedDomain, clonedParticle);
    }
}

template <typename T, template <typename U> class Descriptor>
CloneParticles3D<T, Descriptor> *CloneParticles3D<T, Descriptor>::clone() const
{
    return new CloneParticles3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CloneParticles3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;           // Original particle field.
    modified[1] = modif::dynamicVariables;  // Cloned particle field.
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT CloneParticles3D<T, Descriptor>::appliesTo() const
{
    // The data processor acts on envelope too, but extension to the envelope
    //   is done manually in processGenericBlocks.
    return BlockDomain::bulk;
}

/* ******** CopyParticleToVertex3D *********************************** */

template <typename T, template <typename U> class Descriptor>
CopyParticleToVertex3D<T, Descriptor>::CopyParticleToVertex3D(TriangularSurfaceMesh<T> &mesh_) :
    mesh(mesh_)
{ }

template <typename T, template <typename U> class Descriptor>
void CopyParticleToVertex3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField3D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    // Manually extend the domain to the full envelope.
    domain = particleField.getBoundingBox();

    std::vector<Particle3D<T, Descriptor> *> found;
    particleField.findParticles(domain, found);

    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        Particle3D<T, Descriptor> &particle = *found[iParticle];
        Array<T, 3> position(particle.getPosition());
        plint vertexId = particle.getTag();
        mesh.replaceVertex(vertexId, position);
    }
}

template <typename T, template <typename U> class Descriptor>
CopyParticleToVertex3D<T, Descriptor> *CopyParticleToVertex3D<T, Descriptor>::clone() const
{
    return new CopyParticleToVertex3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT CopyParticleToVertex3D<T, Descriptor>::appliesTo() const
{
    // The data processor acts on envelope too, but extension to the envelope
    //   is done manually in processGenericBlocks.
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void CopyParticleToVertex3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;  // Particle field.
}

/* ******** ComputeParticleForce3D *********************************** */

template <typename T, template <typename U> class Descriptor>
ComputeParticleForce3D<T, Descriptor>::ComputeParticleForce3D(
    TriangleBoundary3D<T> const &triangleBoundary_, int flowType_, bool incompressibleModel_) :
    triangleBoundary(triangleBoundary_),
    flowType(flowType_),
    incompressibleModel(incompressibleModel_)
{ }

template <typename T, template <typename U> class Descriptor>
void ComputeParticleForce3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);
    ParticleField3D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    BlockLattice3D<T, Descriptor> &fluid =
        *dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[1]);
    ScalarField3D<int> &voxelMatrix = *dynamic_cast<ScalarField3D<int> *>(blocks[2]);

    // A few comments on the algorithm. To get the values of a macroscopic quantity
    // on a mesh vertex, the quantity is interpolated from the 4 (2D) or 8 (3D)
    // neighbors. The problem: all these neighbors are not necessarily well defined.
    // In the following scheme, which represents a possible situation in an interior
    // flow, the "o" stand for exterior nodes and the "x" for a boundary vertex. The
    // two left-most "o" are undefined, because none of their neighbors are interior
    // nodes.
    /*
               /
        o    o/
             /
            x
             \
        o    o\
               \
    */
    // These two nodes are therefore manually excluded in the code.

    std::vector<Particle3D<T, Descriptor> *> found;
    particleField.findParticles(domain, found);

    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    std::vector<Cell<T, Descriptor> *> cells(8);
    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        Particle3D<T, Descriptor> &particle = *found[iParticle];
        Array<T, 3> position(particle.getPosition());
        plint vertexId = particle.getTag();
        Array<T, 3> n = triangleBoundary.getMesh().computeVertexNormal(vertexId);

        plint numUsableCells = 0;
        plint numTrials = 0;
        while (numUsableCells == 0 && numTrials < 5) {
            linearInterpolationCoefficients(fluid, position, cellPos, weights);
            for (plint iCell = 0; iCell < 8; ++iCell) {
                cells[iCell] = &fluid.get(cellPos[iCell].x, cellPos[iCell].y, cellPos[iCell].z);
                int flag = voxelMatrix.get(cellPos[iCell].x, cellPos[iCell].y, cellPos[iCell].z);
                // In an interior flow for example, we exclude outer non-boundary cells, because
                //   they have undefined flow data.
                if (flag != voxelFlag::invert(flowType)) {
                    ++numUsableCells;
                }
            }
            position -= (T)1.e-1 * n;
            ++numTrials;
        }
        // PLB_ASSERT(numUsableCells>0);
        if (numUsableCells <= 0)
        {  // TODO: This is a quick fix around cases where no usable cell was found around the
           // vertex.
           //       It would be better to find a more robust handling of the geometry.
            std::vector<Array<T, 3> > vectors;
            Array<T, 3> forceOnWall((T)0., (T)0., (T)0.);
            vectors.push_back(forceOnWall);
            particle.setVectors(vectors);

            T pressure = Descriptor<T>::cs2;
            T wss = 0.0;

            std::vector<T> scalars;
            scalars.push_back(pressure);
            scalars.push_back(wss);
            particle.setScalars(scalars);
            continue;
        }

        Cell<T, Descriptor> *cellOnVertex = 0;
        if (numUsableCells == 8) {
            // Use copy constructor in order to initialize dynamics object.
            cellOnVertex = new Cell<T, Descriptor>(*cells[0]);
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                (*cellOnVertex)[iPop] =
                    weights[0] * (*cells[0])[iPop] + weights[1] * (*cells[1])[iPop]
                    + weights[2] * (*cells[2])[iPop] + weights[3] * (*cells[3])[iPop]
                    + weights[4] * (*cells[4])[iPop] + weights[5] * (*cells[5])[iPop]
                    + weights[6] * (*cells[6])[iPop] + weights[7] * (*cells[7])[iPop];
            }
        } else {
            bool firstUsableCell = true;
            T sumWeights = T();
            for (plint iCell = 0; iCell < 8; ++iCell) {
                int flag = voxelMatrix.get(cellPos[iCell].x, cellPos[iCell].y, cellPos[iCell].z);
                // In an interior flow for example, we exclude outer non-boundary cells, because
                //   they have undefined flow data.
                if (flag != voxelFlag::invert(flowType)) {
                    if (firstUsableCell) {
                        // Use copy constructor in order to initialize dynamics object.
                        cellOnVertex = new Cell<T, Descriptor>(*cells[iCell]);
                        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                            (*cellOnVertex)[iPop] = T();
                        }
                        firstUsableCell = false;
                    }
                    sumWeights += weights[iCell];
                    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                        (*cellOnVertex)[iPop] += weights[iCell] * (*cells[iCell])[iPop];
                    }
                }
            }
            T invSumWeights = (T)1 / sumWeights;
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                (*cellOnVertex)[iPop] *= invSumWeights;
            }
        }

        T rhoBar;
        Array<T, 3> j;
        Array<T, SymmetricTensorImpl<T, 3>::n> PiNeq;
        cellOnVertex->getDynamics().computeRhoBarJPiNeq(*cellOnVertex, rhoBar, j, PiNeq);
        T rho = Descriptor<T>::fullRho(rhoBar);

        Array<T, 3> Pi_n;
        SymmetricTensorImpl<T, 3>::matVectMult(PiNeq, n, Pi_n);

        T omega = cellOnVertex->getDynamics().getOmega();
        Array<T, 3> forceOnFluid;
        if (incompressibleModel) {
            forceOnFluid = -(rho - (T)1.) * Descriptor<T>::cs2 * n
                           + (omega / (T)2. - (T)1.) * Pi_n;  // Incompressible vision
        } else {
            forceOnFluid = -(rho - (T)1.) * Descriptor<T>::cs2 * n
                           + Descriptor<T>::invRho(rhoBar) * (omega / (T)2. - (T)1.)
                                 * Pi_n;  // Compressible vision
        }
        Array<T, 3> forceOnWall = -forceOnFluid;  // Remark by J.L: This is Newton!
        delete cellOnVertex;

        std::vector<Array<T, 3> > vectors;
        vectors.push_back(forceOnWall);
        particle.setVectors(vectors);

        T pressure = (rho - (T)1.) * Descriptor<T>::cs2;
        T wss = norm(forceOnWall - dot(n, forceOnWall) * n);

        std::vector<T> scalars;
        scalars.push_back(pressure);
        scalars.push_back(wss);
        particle.setScalars(scalars);
    }
}

template <typename T, template <typename U> class Descriptor>
ComputeParticleForce3D<T, Descriptor> *ComputeParticleForce3D<T, Descriptor>::clone() const
{
    return new ComputeParticleForce3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ComputeParticleForce3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;           // Fluid.
    modified[2] = modif::nothing;           // Voxel Matrix.
}

/* ******** GetParticleVelocity3D *********************************** */

template <typename T, template <typename U> class Descriptor>
void GetParticleVelocity3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField3D<T, Descriptor> &originalField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    ParticleField3D<T, Descriptor> &clonedField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[1]);

    std::vector<Particle3D<T, Descriptor> *> found;
    originalField.findParticles(domain, found);

    Dot3D offset = computeRelativeDisplacement(originalField, clonedField);
    Box3D clonedDomain(domain.shift(offset.x, offset.y, offset.z));

    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        Particle3D<T, Descriptor> &originalParticle = *found[iParticle];
        Particle3D<T, Descriptor> *clonedParticle = new VisualParticle3D<T, Descriptor>(
            originalParticle.getTag(), originalParticle.getPosition());
        std::vector<Array<T, 3> > vectors;
        Array<T, 3> velocity;
#ifdef PLB_DEBUG
        bool ok =
#endif
            originalParticle.getVector(0, velocity);
        PLB_ASSERT(ok);
        vectors.push_back(velocity);
        clonedParticle->setVectors(vectors), clonedField.addParticle(clonedDomain, clonedParticle);
    }
}

template <typename T, template <typename U> class Descriptor>
GetParticleVelocity3D<T, Descriptor> *GetParticleVelocity3D<T, Descriptor>::clone() const
{
    return new GetParticleVelocity3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT GetParticleVelocity3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;  // Important: must be bulk-only,
    // because non-local access to fluid is made.
}

template <typename T, template <typename U> class Descriptor>
void GetParticleVelocity3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;           // Original field.
    modified[1] = modif::dynamicVariables;  // Cloned field.
}

/* ******** CreateParticleMap3D *********************************** */

template <typename T, template <typename U> class Descriptor>
void CreateParticleMap3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 4);
    ParticleField3D<T, Descriptor> &particleField =
        dynamic_cast<ParticleField3D<T, Descriptor> &>(*blocks[0]);
    ScalarField3D<T> &xMatrix = dynamic_cast<ScalarField3D<T> &>(*blocks[1]);
    ScalarField3D<T> &yMatrix = dynamic_cast<ScalarField3D<T> &>(*blocks[2]);
    ScalarField3D<T> &zMatrix = dynamic_cast<ScalarField3D<T> &>(*blocks[3]);

    ScalarField3D<int> counter(xMatrix.getNx(), xMatrix.getNy(), xMatrix.getNz(), 0);

    std::vector<Particle3D<T, Descriptor> *> found;
    particleField.findParticles(domain, found);

    Dot3D location(particleField.getLocation());
    Dot3D offset = computeRelativeDisplacement(particleField, xMatrix);

    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        Particle3D<T, Descriptor> *originalParticle = found[iParticle];
        MappingParticle3D<T, Descriptor> *mappingParticle =
            dynamic_cast<MappingParticle3D<T, Descriptor> *>(originalParticle);
        if (mappingParticle && mappingParticle->endOfCycle()) {
            Array<T, 3> pos1(mappingParticle->getPos1());
            Array<T, 3> pos2(mappingParticle->getPos2());
            plint xPos = (plint)(pos1[0]) - location.x;
            plint yPos = (plint)(pos1[1]) - location.y;
            plint zPos = (plint)(pos1[2]) - location.z;
            plint xPosOfs = xPos + offset.x;
            plint yPosOfs = yPos + offset.y;
            plint zPosOfs = zPos + offset.z;
            ++counter.get(xPosOfs, yPosOfs, zPosOfs);
            xMatrix.get(xPosOfs, yPosOfs, zPosOfs) += pos2[0];
            yMatrix.get(xPosOfs, yPosOfs, zPosOfs) += pos2[1];
            zMatrix.get(xPosOfs, yPosOfs, zPosOfs) += pos2[2];
        }
    }
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint xPosOfs = iX + offset.x;
                plint yPosOfs = iY + offset.y;
                plint zPosOfs = iZ + offset.z;
                int count = counter.get(xPosOfs, yPosOfs, zPosOfs);
                if (count == 0) {
                    xMatrix.get(xPosOfs, yPosOfs, zPosOfs) = -1.;
                    yMatrix.get(xPosOfs, yPosOfs, zPosOfs) = -1.;
                    zMatrix.get(xPosOfs, yPosOfs, zPosOfs) = -1.;
                } else {
                    xMatrix.get(xPosOfs, yPosOfs, zPosOfs) /= (T)count;
                    yMatrix.get(xPosOfs, yPosOfs, zPosOfs) /= (T)count;
                    zMatrix.get(xPosOfs, yPosOfs, zPosOfs) /= (T)count;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CreateParticleMap3D<T, Descriptor> *CreateParticleMap3D<T, Descriptor>::clone() const
{
    return new CreateParticleMap3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CreateParticleMap3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Particle field.
    modified[1] = modif::staticVariables;  // x-Matrix.
    modified[2] = modif::staticVariables;  // y-Matrix.
    modified[3] = modif::staticVariables;  // z-Matrix.
}

/* ******** VisualScalarFieldParticles3D *********************************** */

template <typename T, template <typename U> class Descriptor>
VisualScalarFieldParticles3D<T, Descriptor>::VisualScalarFieldParticles3D(
    T mostUnlikely_, T mostLikely_, plint numShotsPerCell_, T probability_, plint ny_, plint nz_,
    Particle3D<T, Descriptor> *particleTemplate_) :
    mostUnlikely(mostUnlikely_),
    mostLikely(mostLikely_),
    numShotsPerCell(numShotsPerCell_),
    probability(probability_),
    ny(ny_),
    nz(nz_),
    particleTemplate(particleTemplate_)
{ }

template <typename T, template <typename U> class Descriptor>
VisualScalarFieldParticles3D<T, Descriptor>::~VisualScalarFieldParticles3D()
{
    delete particleTemplate;
}

template <typename T, template <typename U> class Descriptor>
VisualScalarFieldParticles3D<T, Descriptor>::VisualScalarFieldParticles3D(
    VisualScalarFieldParticles3D<T, Descriptor> const &rhs) :
    mostUnlikely(rhs.mostUnlikely),
    mostLikely(rhs.mostLikely),
    numShotsPerCell(rhs.numShotsPerCell),
    probability(rhs.probability),
    ny(rhs.ny),
    nz(rhs.nz),
    particleTemplate(rhs.particleTemplate->clone())

{ }

template <typename T, template <typename U> class Descriptor>
VisualScalarFieldParticles3D<T, Descriptor> &VisualScalarFieldParticles3D<T, Descriptor>::operator=(
    VisualScalarFieldParticles3D<T, Descriptor> const &rhs)
{
    VisualScalarFieldParticles3D<T, Descriptor>(rhs).swap(*this);
}

template <typename T, template <typename U> class Descriptor>
void VisualScalarFieldParticles3D<T, Descriptor>::swap(
    VisualScalarFieldParticles3D<T, Descriptor> &rhs)
{
    std::swap(mostUnlikely, rhs.mostUnlikely);
    std::swap(mostLikely, rhs.mostLikely);
    std::swap(numShotsPerCell, rhs.numShotsPerCell);
    std::swap(probability, rhs.probability);
    std::swap(ny, rhs.ny);
    std::swap(nz, rhs.nz);
    std::swap(particleTemplate, rhs.particleTemplate);
}

template <typename T, template <typename U> class Descriptor>
void VisualScalarFieldParticles3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField3D<T, Descriptor> *particleFieldPtr =
        dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(particleFieldPtr);
    ParticleField3D<T, Descriptor> &particleField = *particleFieldPtr;

    ScalarField3D<T> *scalarPtr = dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    PLB_ASSERT(scalarPtr);
    ScalarField3D<T> &scalarField = *scalarPtr;

    Dot3D const &location = particleField.getLocation();

    sitmo::prng_engine eng;
    plint rng_index = 0;

    T invDelta = mostLikely - mostUnlikely;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint xOffset = ny * (iX + location.x);
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint yOffset = nz * (xOffset + iY + location.y);
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint zOffset = yOffset + iZ + location.z;

                PLB_ASSERT(zOffset >= rng_index);
                if (zOffset > rng_index) {
                    eng.discard(zOffset - rng_index);
                    rng_index = zOffset;
                }
                T randomValue = (T)eng() / ((T)sitmo::prng_engine::max() + 1.0);
                ++rng_index;

                sitmo::prng_engine eng2(zOffset);
                for (plint iTimes = 0; iTimes < numShotsPerCell; ++iTimes) {
                    T posX =
                        (T)iX
                        + (T)0.5 * (randomValue + (T)eng2() / ((T)sitmo::prng_engine::max() + 1.0))
                        - (T)0.5 + location.x;
                    T posY =
                        (T)iY
                        + (T)0.5 * (randomValue + (T)eng2() / ((T)sitmo::prng_engine::max() + 1.0))
                        - (T)0.5 + location.y;
                    T posZ =
                        (T)iZ
                        + (T)0.5 * (randomValue + (T)eng2() / ((T)sitmo::prng_engine::max() + 1.0))
                        - (T)0.5 + location.z;
                    Array<T, 3> pos(posX, posY, posZ);

                    T scalarValue = linearInterpolateScalarField(scalarField, pos);
                    T scaledScalarValue = (scalarValue - mostUnlikely) * invDelta;
                    if (scaledScalarValue < T())
                        scaledScalarValue = T();
                    if (scaledScalarValue > (T)1.)
                        scaledScalarValue = (T)1.;

                    T overallProbability = probability * scaledScalarValue;
                    if ((T)0.5 * (randomValue + (T)eng2() / ((T)sitmo::prng_engine::max() + 1.0))
                        > overallProbability) {
                        Particle3D<T, Descriptor> *newParticle = particleTemplate->clone();
                        newParticle->getPosition() = pos;
                        std::vector<T> scalarVector;
                        scalarVector.push_back(scalarValue);
                        newParticle->setScalars(scalarVector);
                        particleField.addParticle(domain, newParticle);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
VisualScalarFieldParticles3D<T, Descriptor> *VisualScalarFieldParticles3D<T, Descriptor>::clone()
    const
{
    return new VisualScalarFieldParticles3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void VisualScalarFieldParticles3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;           // Scalar field.
}

/* ******** ResetInitialVelocityOfPlaneMappingParticles3D *********************************** */

template <typename T, template <typename U> class Descriptor>
void ResetInitialVelocityOfPlaneMappingParticles3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 1);
    ParticleField3D<T, Descriptor> *particleField =
        dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(particleField);

    std::vector<Particle3D<T, Descriptor> *> particles;
    particleField->findParticles(domain, particles);
    for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
        PlaneMappingParticle3D<T, Descriptor> *particle =
            dynamic_cast<PlaneMappingParticle3D<T, Descriptor> *>(particles[iParticle]);
        PLB_ASSERT(particle);
        particle->getInitialVelocity() = particle->getVelocity();
    }
}

template <typename T, template <typename U> class Descriptor>
ResetInitialVelocityOfPlaneMappingParticles3D<T, Descriptor>
    *ResetInitialVelocityOfPlaneMappingParticles3D<T, Descriptor>::clone() const
{
    return new ResetInitialVelocityOfPlaneMappingParticles3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ResetInitialVelocityOfPlaneMappingParticles3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field of PlaneMappingParticle3D particles.
}

/* ******** AssignTimeScalingToPlaneMappingParticles3D *********************************** */

template <typename T, template <typename U> class Descriptor>
AssignTimeScalingToPlaneMappingParticles3D<
    T, Descriptor>::AssignTimeScalingToPlaneMappingParticles3D(T timeScaling_) :
    timeScaling(timeScaling_)
{ }

template <typename T, template <typename U> class Descriptor>
void AssignTimeScalingToPlaneMappingParticles3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 1);
    ParticleField3D<T, Descriptor> *particleField =
        dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(particleField);

    std::vector<Particle3D<T, Descriptor> *> particles;
    particleField->findParticles(domain, particles);
    for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
        PlaneMappingParticle3D<T, Descriptor> *particle =
            dynamic_cast<PlaneMappingParticle3D<T, Descriptor> *>(particles[iParticle]);
        PLB_ASSERT(particle);
        particle->getTimeScaling() = timeScaling;
    }
}

template <typename T, template <typename U> class Descriptor>
AssignTimeScalingToPlaneMappingParticles3D<T, Descriptor>
    *AssignTimeScalingToPlaneMappingParticles3D<T, Descriptor>::clone() const
{
    return new AssignTimeScalingToPlaneMappingParticles3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AssignTimeScalingToPlaneMappingParticles3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field of PlaneMappingParticle3D particles.
}

/* ******** PlaneMappingParticleStatistics3D *********************************** */

template <typename T, template <typename U> class Descriptor>
PlaneMappingParticleStatistics3D<T, Descriptor>::PlaneMappingParticleStatistics3D(
    T cutOffValue_, int wallFlag_, int fluidFlag_) :
    cutOffValue(cutOffValue_),
    wallFlag(wallFlag_),
    fluidFlag(fluidFlag_),
    numParticlesId(this->getStatistics().subscribeIntSum()),
    numTerminalParticlesId(this->getStatistics().subscribeIntSum()),
    numWallParticlesId(this->getStatistics().subscribeIntSum()),
    numSlowParticlesId(this->getStatistics().subscribeIntSum()),
    numSlowNearWallParticlesId(this->getStatistics().subscribeIntSum())
{ }

template <typename T, template <typename U> class Descriptor>
void PlaneMappingParticleStatistics3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 2);
    ParticleField3D<T, Descriptor> *particleField =
        dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(particleField);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int> *>(blocks[1]);
    PLB_ASSERT(flag);

    Dot3D ofs = computeRelativeDisplacement(*particleField, *flag);

    plint numParticles = 0;
    plint numTerminalParticles = 0;
    plint numWallParticles = 0;
    plint numSlowParticles = 0;
    plint numSlowNearWallParticles = 0;

    std::vector<Particle3D<T, Descriptor> *> particles;
    particleField->findParticles(domain, particles);
    if (particles.size() == 0) {
        return;
    }
    numParticles = (plint)particles.size();
    for (pluint iParticle = 0; iParticle < (pluint)particles.size(); iParticle++) {
        PlaneMappingParticle3D<T, Descriptor> *particle =
            dynamic_cast<PlaneMappingParticle3D<T, Descriptor> *>(particles[iParticle]);
        PLB_ASSERT(particle);
        if (particle->passedTerminalPlane()) {
            numTerminalParticles++;
        }
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (flag->get(iX + ofs.x, iY + ofs.y, iZ + ofs.z) == wallFlag) {
                    Box3D node(iX, iX, iY, iY, iZ, iZ);
                    std::vector<Particle3D<T, Descriptor> *> wallParticles;
                    particleField->findParticles(node, wallParticles);
                    if (wallParticles.size() == 0) {
                        continue;
                    }

                    int closeToFluid = 0;
                    for (int i = -1; i < 2; i++) {
                        for (int j = -1; j < 2; j++) {
                            for (int k = -1; k < 2; k++) {
                                if (flag->get(iX + i + ofs.x, iY + j + ofs.y, iZ + k + ofs.z)
                                    == fluidFlag) {
                                    closeToFluid = 1;
                                    goto closeToFluidLabel;
                                }
                            }
                        }
                    }
closeToFluidLabel:
                    if (!closeToFluid) {
                        numWallParticles += (plint)wallParticles.size();
                    }
                } else if (flag->get(iX + ofs.x, iY + ofs.y, iZ + ofs.z) == fluidFlag) {
                    Box3D node(iX, iX, iY, iY, iZ, iZ);
                    std::vector<Particle3D<T, Descriptor> *> fluidParticles;
                    particleField->findParticles(node, fluidParticles);
                    if (fluidParticles.size() == 0) {
                        continue;
                    }
                    plint slowParticles = 0;
                    for (pluint iParticle = 0; iParticle < (pluint)fluidParticles.size();
                         iParticle++) {
                        PlaneMappingParticle3D<T, Descriptor> *particle =
                            dynamic_cast<PlaneMappingParticle3D<T, Descriptor> *>(
                                particles[iParticle]);
                        PLB_ASSERT(particle);
                        if (normSqr(particle->getVelocity()) <= cutOffValue) {
                            slowParticles++;
                        }
                    }
                    numSlowParticles += slowParticles;

                    int closeToWall = 0;
                    for (int i = -1; i < 2; i++) {
                        for (int j = -1; j < 2; j++) {
                            for (int k = -1; k < 2; k++) {
                                if (flag->get(iX + i + ofs.x, iY + j + ofs.y, iZ + k + ofs.z)
                                    == wallFlag) {
                                    closeToWall = 1;
                                    goto closeToWallLabel;
                                }
                            }
                        }
                    }
closeToWallLabel:
                    if (closeToWall) {
                        numSlowNearWallParticles += slowParticles;
                    }
                }
            }
        }
    }

    this->getStatistics().gatherIntSum(numParticlesId, numParticles);
    this->getStatistics().gatherIntSum(numTerminalParticlesId, numTerminalParticles);
    this->getStatistics().gatherIntSum(numWallParticlesId, numWallParticles);
    this->getStatistics().gatherIntSum(numSlowParticlesId, numSlowParticles);
    this->getStatistics().gatherIntSum(numSlowNearWallParticlesId, numSlowNearWallParticles);
}

template <typename T, template <typename U> class Descriptor>
PlaneMappingParticleStatistics3D<T, Descriptor>
    *PlaneMappingParticleStatistics3D<T, Descriptor>::clone() const
{
    return new PlaneMappingParticleStatistics3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void PlaneMappingParticleStatistics3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;  // Particle field.
    modified[1] = modif::nothing;  // Flag matrix.
}

template <typename T, template <typename U> class Descriptor>
plint PlaneMappingParticleStatistics3D<T, Descriptor>::getNumParticles() const
{
    return this->getStatistics().getIntSum(numParticlesId);
}

template <typename T, template <typename U> class Descriptor>
plint PlaneMappingParticleStatistics3D<T, Descriptor>::getNumTerminalParticles() const
{
    return this->getStatistics().getIntSum(numTerminalParticlesId);
}

template <typename T, template <typename U> class Descriptor>
plint PlaneMappingParticleStatistics3D<T, Descriptor>::getNumWallParticles() const
{
    return this->getStatistics().getIntSum(numWallParticlesId);
}

template <typename T, template <typename U> class Descriptor>
plint PlaneMappingParticleStatistics3D<T, Descriptor>::getNumSlowParticles() const
{
    return this->getStatistics().getIntSum(numSlowParticlesId);
}

template <typename T, template <typename U> class Descriptor>
plint PlaneMappingParticleStatistics3D<T, Descriptor>::getNumSlowNearWallParticles() const
{
    return this->getStatistics().getIntSum(numSlowNearWallParticlesId);
}

template <typename T, template <typename U> class Descriptor>
void computePlaneMappingParticleStatistics3D(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles,
    MultiScalarField3D<int> &flagMatrix, T cutOffValue, int wallFlag, int fluidFlag, Box3D domain,
    plint &numParticles, plint &numTerminalParticles, plint &numWallParticles,
    plint &numSlowParticles, plint &numSlowNearWallParticles)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&particles);
    args.push_back(&flagMatrix);
    PlaneMappingParticleStatistics3D<T, Descriptor> functional(cutOffValue, wallFlag, fluidFlag);
    applyProcessingFunctional(functional, domain, args);

    numParticles = functional.getNumParticles();
    numTerminalParticles = functional.getNumTerminalParticles();
    numWallParticles = functional.getNumWallParticles();
    numSlowParticles = functional.getNumSlowParticles();
    numSlowNearWallParticles = functional.getNumSlowNearWallParticles();
}

/* ******** LinearScalarFieldToParticle3D *********************************** */

template <typename T, template <typename U> class Descriptor>
LinearScalarFieldToParticle3D<T, Descriptor>::LinearScalarFieldToParticle3D(
    int flowType_, plint whichScalar_) :
    whichScalar(whichScalar_)
{
    usableFlowTypes.push_back(flowType_);
}

template <typename T, template <typename U> class Descriptor>
LinearScalarFieldToParticle3D<T, Descriptor>::LinearScalarFieldToParticle3D(
    std::vector<int> const &usableFlowTypes_, plint whichScalar_) :
    usableFlowTypes(usableFlowTypes_), whichScalar(whichScalar_)
{ }

template <typename T, template <typename U> class Descriptor>
void LinearScalarFieldToParticle3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);
    ParticleField3D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    ScalarField3D<T> &scalarField = *dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    ScalarField3D<int> &voxelMatrix = *dynamic_cast<ScalarField3D<int> *>(blocks[2]);

    Dot3D location = particleField.getLocation();
    Dot3D ofs = computeRelativeDisplacement(particleField, scalarField);
    Dot3D ofsV = computeRelativeDisplacement(particleField, voxelMatrix);

    std::vector<Particle3D<T, Descriptor> *> found;
    particleField.findParticles(domain, found);

    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        Particle3D<T, Descriptor> &particle = *found[iParticle];
        Array<T, 3> vertex = particle.getPosition();
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
        T averageScalar = T();
        T totWeight = T();
        // x   x . x   x
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                    T nextScalar = scalarField.get(pos[0] + ofs.x, pos[1] + ofs.y, pos[2] + ofs.z);

                    T nextFlag = voxelMatrix.get(pos[0] + ofsV.x, pos[1] + ofsV.y, pos[2] + ofsV.z);

                    std::vector<int>::iterator it =
                        std::find(usableFlowTypes.begin(), usableFlowTypes.end(), nextFlag);
                    if (it != usableFlowTypes.end()) {
                        Array<T, 3> r((T)dx - fracPos[0], (T)dy - fracPos[1], (T)dz - fracPos[2]);
                        T W = inamuroDeltaFunction<T>().W(r);
                        averageScalar += W * nextScalar;
                        totWeight += W;
                    }
                }
            }
        }
        averageScalar /= totWeight;

        if (whichScalar == -1) {
            std::vector<T> scalars;
            scalars.push_back(averageScalar);
            particle.setScalars(scalars);
        } else {
            particle.setScalar(whichScalar, averageScalar);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LinearScalarFieldToParticle3D<T, Descriptor> *LinearScalarFieldToParticle3D<T, Descriptor>::clone()
    const
{
    return new LinearScalarFieldToParticle3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void LinearScalarFieldToParticle3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;           // Scalar field.
    modified[2] = modif::nothing;           // Scalar field.
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT LinearScalarFieldToParticle3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ScalarFieldToParticle3D *********************************** */

template <typename T, template <typename U> class Descriptor>
ScalarFieldToParticle3D<T, Descriptor>::ScalarFieldToParticle3D(plint whichScalar_) :
    whichScalar(whichScalar_)
{ }

template <typename T, template <typename U> class Descriptor>
void ScalarFieldToParticle3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField3D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    ScalarField3D<T> &scalarField = *dynamic_cast<ScalarField3D<T> *>(blocks[1]);

    Dot3D location = particleField.getLocation();
    Dot3D ofs = computeRelativeDisplacement(particleField, scalarField);

    std::vector<Particle3D<T, Descriptor> *> found;
    particleField.findParticles(domain, found);

    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        Particle3D<T, Descriptor> &particle = *found[iParticle];
        Array<T, 3> vertex = particle.getPosition();
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
        T averageScalar = T();
        // x   x . x   x
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                    T nextScalar = scalarField.get(pos[0] + ofs.x, pos[1] + ofs.y, pos[2] + ofs.z);
                    Array<T, 3> r((T)dx - fracPos[0], (T)dy - fracPos[1], (T)dz - fracPos[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageScalar += W * nextScalar;
                }
            }
        }

        if (whichScalar == -1) {
            std::vector<T> scalars;
            scalars.push_back(averageScalar);
            particle.setScalars(scalars);
        } else {
            particle.setScalar(whichScalar, averageScalar);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ScalarFieldToParticle3D<T, Descriptor> *ScalarFieldToParticle3D<T, Descriptor>::clone() const
{
    return new ScalarFieldToParticle3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ScalarFieldToParticle3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;           // Scalar field.
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT ScalarFieldToParticle3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** TensorFieldToParticle3D *********************************** */

template <typename T, int nDim, template <typename U> class Descriptor>
TensorFieldToParticle3D<T, nDim, Descriptor>::TensorFieldToParticle3D(plint whichScalar_) :
    whichScalar(whichScalar_)
{ }

template <typename T, int nDim, template <typename U> class Descriptor>
void TensorFieldToParticle3D<T, nDim, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField3D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    TensorField3D<T, nDim> &tensorField = *dynamic_cast<TensorField3D<T, nDim> *>(blocks[1]);

    Dot3D location = particleField.getLocation();
    Dot3D ofs = computeRelativeDisplacement(particleField, tensorField);

    std::vector<Particle3D<T, Descriptor> *> found;
    particleField.findParticles(domain, found);

    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        Particle3D<T, Descriptor> &particle = *found[iParticle];
        Array<T, 3> vertex = particle.getPosition();
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
        std::vector<T> averageData(nDim, T());
        // x   x . x   x
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                    Array<T, nDim> nextData =
                        tensorField.get(pos[0] + ofs.x, pos[1] + ofs.y, pos[2] + ofs.z);
                    Array<T, 3> r((T)dx - fracPos[0], (T)dy - fracPos[1], (T)dz - fracPos[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    for (plint i = 0; i < nDim; ++i) {
                        averageData[i] += W * nextData[i];
                    }
                }
            }
        }

        if (whichScalar == -1) {
            particle.setScalars(averageData);
        } else {
            for (plint i = 0; i < nDim; ++i) {
                particle.setScalar(whichScalar + i, averageData[i]);
            }
        }
    }
}

template <typename T, int nDim, template <typename U> class Descriptor>
TensorFieldToParticle3D<T, nDim, Descriptor> *TensorFieldToParticle3D<T, nDim, Descriptor>::clone()
    const
{
    return new TensorFieldToParticle3D<T, nDim, Descriptor>(*this);
}

template <typename T, int nDim, template <typename U> class Descriptor>
void TensorFieldToParticle3D<T, nDim, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;           // Tensor field.
}

}  // namespace plb

#endif  // VISUAL_PARTICLE_FUNCTIONAL_3D_HH
