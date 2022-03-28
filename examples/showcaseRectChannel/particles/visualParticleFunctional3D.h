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

#ifndef VISUAL_PARTICLE_FUNCTIONAL_3D_H
#define VISUAL_PARTICLE_FUNCTIONAL_3D_H

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/globalDefs.h"
#include "offLattice/triangleBoundary3D.h"
#include "offLattice/triangularSurfaceMesh.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, class ParticleT>
class CreateParticleFromVertex3D : public BoxProcessingFunctional3D {
public:
    CreateParticleFromVertex3D(TriangularSurfaceMesh<T> const &mesh_);
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CreateParticleFromVertex3D<T, Descriptor, ParticleT> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    TriangularSurfaceMesh<T> const &mesh;
};

template <typename T, template <typename U> class Descriptor>
class CreateParticleFromVertexNoTagging3D : public BoxProcessingFunctional3D {
public:
    CreateParticleFromVertexNoTagging3D(
        TriangularSurfaceMesh<T> const &mesh_, Particle3D<T, Descriptor> *particleTemplate_);
    ~CreateParticleFromVertexNoTagging3D();
    CreateParticleFromVertexNoTagging3D(
        CreateParticleFromVertexNoTagging3D<T, Descriptor> const &rhs);
    CreateParticleFromVertexNoTagging3D<T, Descriptor> &operator=(
        CreateParticleFromVertexNoTagging3D<T, Descriptor> const &rhs);
    void swap(CreateParticleFromVertexNoTagging3D<T, Descriptor> &rhs);
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CreateParticleFromVertexNoTagging3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    TriangularSurfaceMesh<T> const *mesh;
    Particle3D<T, Descriptor> *particleTemplate;
};

template <typename T, template <typename U> class Descriptor>
class CloneParticles3D : public BoxProcessingFunctional3D {
public:
    CloneParticles3D(Particle3D<T, Descriptor> *particleTemplate_);
    CloneParticles3D(CloneParticles3D<T, Descriptor> const &rhs);
    ~CloneParticles3D();
    CloneParticles3D<T, Descriptor> &operator=(CloneParticles3D<T, Descriptor> const &rhs);
    void swap(CloneParticles3D<T, Descriptor> &rhs);
    /// Arguments: [0] original particles; [1] cloned particles.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CloneParticles3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
};

template <typename T, template <typename U> class Descriptor, class ParticleT>
class CopyVertexToParticle3D : public BoxProcessingFunctional3D {
public:
    CopyVertexToParticle3D(TriangularSurfaceMesh<T> const &mesh_);
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CopyVertexToParticle3D<T, Descriptor, ParticleT> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    TriangularSurfaceMesh<T> const &mesh;
};

template <typename T, template <typename U> class Descriptor>
class CopyParticleToVertex3D : public BoxProcessingFunctional3D {
public:
    CopyParticleToVertex3D(TriangularSurfaceMesh<T> &mesh_);
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CopyParticleToVertex3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    TriangularSurfaceMesh<T> &mesh;
};

template <typename T, template <typename U> class Descriptor>
class ComputeParticleForce3D : public BoxProcessingFunctional3D {
public:
    ComputeParticleForce3D(
        TriangleBoundary3D<T> const &triangleBoundary_, int flowType_, bool incompressibleModel_);
    /// Arguments: [0] Particle-field; [1] Scalar-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ComputeParticleForce3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    TriangleBoundary3D<T> const &triangleBoundary;
    int flowType;
    bool incompressibleModel;
};

template <typename T, template <typename U> class Descriptor>
class LinearScalarFieldToParticle3D : public BoxProcessingFunctional3D {
public:
    LinearScalarFieldToParticle3D(int flowType_, plint whichScalar_ = -1);
    LinearScalarFieldToParticle3D(
        std::vector<int> const &usableFlowTypes_, plint whichScalar_ = -1);
    /// Arguments: [0] Particle-field; [1] Scalar-Field; [2] Voxel matrix
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual LinearScalarFieldToParticle3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    std::vector<int> usableFlowTypes;
    plint whichScalar;
};

template <typename T, template <typename U> class Descriptor>
class ScalarFieldToParticle3D : public BoxProcessingFunctional3D {
public:
    ScalarFieldToParticle3D(plint whichScalar_ = -1);
    /// Arguments: [0] Particle-field; [1] Scalar-Field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ScalarFieldToParticle3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    plint whichScalar;
};

template <typename T, int nDim, template <typename U> class Descriptor>
class TensorFieldToParticle3D : public BoxProcessingFunctional3D {
public:
    TensorFieldToParticle3D(plint whichScalar_ = -1);
    /// Arguments: [0] Particle-field; [1] Scalar-Field
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual TensorFieldToParticle3D<T, nDim, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint whichScalar;
};

template <typename T, template <typename U> class Descriptor>
class GetParticleVelocity3D : public BoxProcessingFunctional3D {
public:
    /// Arguments: [0] Orignal field (point-particles), [1] Cloned field (visual particles).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual GetParticleVelocity3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, template <typename U> class Descriptor>
class CreateParticleMap3D : public BoxProcessingFunctional3D {
public:
    /// Arguments: [0] particle field; [1] x-matrix of map; [2] y-matrix of map; [3] z-matrix of
    /// map.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CreateParticleMap3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/// Create particles at weighted random positions for visualization purposes.
template <typename T, template <typename U> class Descriptor>
class VisualScalarFieldParticles3D : public BoxProcessingFunctional3D {
public:
    VisualScalarFieldParticles3D(
        T mostUnlikely_, T mostLikely_, plint numShotsPerCell_, T probability_, plint ny_,
        plint nz_, Particle3D<T, Descriptor> *particleTemplate_);
    ~VisualScalarFieldParticles3D();
    VisualScalarFieldParticles3D(VisualScalarFieldParticles3D<T, Descriptor> const &rhs);
    VisualScalarFieldParticles3D<T, Descriptor> &operator=(
        VisualScalarFieldParticles3D<T, Descriptor> const &rhs);
    void swap(VisualScalarFieldParticles3D<T, Descriptor> &rhs);
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual VisualScalarFieldParticles3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T mostUnlikely, mostLikely;
    plint numShotsPerCell;
    T probability;
    plint ny, nz;
    Particle3D<T, Descriptor> *particleTemplate;
};

/// Take a particle field of PlaneMappingParticle3D particles, and reset the initial velocity member
/// of each particle to be its current velocity.
template <typename T, template <typename U> class Descriptor>
class ResetInitialVelocityOfPlaneMappingParticles3D : public BoxProcessingFunctional3D {
public:
    /// Arguments: [0] Particle-field of PlaneMappingParticle3D particles.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ResetInitialVelocityOfPlaneMappingParticles3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/// Take a particle field of PlaneMappingParticle3D particles, and assign a time scaling to each one
/// of them.
template <typename T, template <typename U> class Descriptor>
class AssignTimeScalingToPlaneMappingParticles3D : public BoxProcessingFunctional3D {
public:
    AssignTimeScalingToPlaneMappingParticles3D(T timeScaling_);
    /// Arguments: [0] Particle-field of PlaneMappingParticle3D particles.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual AssignTimeScalingToPlaneMappingParticles3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T timeScaling;
};

/// Compute statistics for PlaneMappingParticles3D
template <typename T, template <typename U> class Descriptor>
class PlaneMappingParticleStatistics3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    PlaneMappingParticleStatistics3D(T cutOffValue_, int wallFlag_, int fluidFlag_);
    /// Arguments: [0] Particle-field [1] Flag-matrix
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual PlaneMappingParticleStatistics3D<T, Descriptor> *clone() const;
    plint getNumParticles() const;          // Total number of particles.
    plint getNumTerminalParticles() const;  // Number of particles that reached the terminal plane.
    plint getNumWallParticles() const;      // Number of particles inside walls (if any).
    plint getNumSlowParticles() const;      // Number of particles with very small velocities.
    plint getNumSlowNearWallParticles() const;  // Number of particles with very small velocities
                                                // which are in the viscinity of walls.
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T cutOffValue;  // When the speed of the particle drops below sqrt(cutOffValue), then this
                    // particle is a "slow" particle.
    int wallFlag;   // Value that represents the wall nodes in the flag matrix.
    int fluidFlag;  // Value that represents the fluid nodes in the flag matrix.
    plint numParticlesId;
    plint numTerminalParticlesId;
    plint numWallParticlesId;
    plint numSlowParticlesId;
    plint numSlowNearWallParticlesId;
};

template <typename T, template <typename U> class Descriptor>
void computePlaneMappingParticleStatistics3D(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles,
    MultiScalarField3D<int> &flagMatrix, T cutOffValue, int wallFlag, int fluidFlag, Box3D domain,
    plint &numParticles, plint &numTerminalParticles, plint &numWallParticles,
    plint &numSlowParticles, plint &numSlowNearWallParticles);

}  // namespace plb

#endif  // VISUAL_PARTICLE_FUNCTIONAL_3D_H
