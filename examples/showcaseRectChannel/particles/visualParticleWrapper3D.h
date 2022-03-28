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

#ifndef VISUAL_PARTICLE_WRAPPER_3D_H
#define VISUAL_PARTICLE_WRAPPER_3D_H

#include <memory>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "offLattice/offLatticeBoundaryCondition3D.h"
#include "particles/visualParticleFunctional3D.h"

namespace plb {

/// Create particles that contain the surface force.
/** Remember: a mesh is a non-parallel object, whereas particles are parallel
 *  objects. This function does not have parallel efficiency, because every
 *  processor needs to loop over the full mesh in order to create the parallel
 *  particles. If you already have the particles, use the version of
 *  computeSurfaceForce() that takes them as an argument, and you will get
 *  parallel efficiency.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > > computeSurfaceForce(
    TriangleBoundary3D<T> &boundary, VoxelizedDomain3D<T> &voxelizedDomain,
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, bool incompressibleModel,
    bool dynamicMesh = false);

/// Create particles that contain the surface force.
/** This version acts on the full domain of the lattice.  **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > > computeSurfaceForce(
    TriangleBoundary3D<T> &boundary, VoxelizedDomain3D<T> &voxelizedDomain,
    MultiBlockLattice3D<T, Descriptor> &lattice, bool incompressibleModel,
    bool dynamicMesh = false);

template <typename T, template <typename U> class Descriptor>
/// Create particles that contain the surface force.
/** This version supposes that you already possess a parallel representation
 *  of the wall in the form of particles. It's faster than the previous
 *  version of computeSurfaceForce().
 **/
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > > computeSurfaceForce(
    TriangleBoundary3D<T> &boundary,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &parallelWall,
    VoxelizedDomain3D<T> &voxelizedDomain, MultiBlockLattice3D<T, Descriptor> &lattice,
    bool incompressibleModel, bool dynamicMesh = false);

/// Create particles that contain the surface force.
/** Remember: a mesh is a non-parallel object, whereas particles are parallel
 *  objects. This function does not have parallel efficiency, because every
 *  processor needs to loop over the full mesh in order to create the parallel
 *  particles. If you already have the particles, use the version of
 *  computeSurfaceForce() that takes them as an argument, and you will get
 *  parallel efficiency.
 **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > > computeSurfaceScalar(
    TriangleBoundary3D<T> &boundary, VoxelizedDomain3D<T> &voxelizedDomain,
    MultiScalarField3D<T> &scalarField, Box3D domain, const std::vector<int> &usableFlowTypes,
    plint whichScalar = -1, bool dynamicMesh = false);

/// Create particles that contain the surface force.
/** This version acts on the full domain of the lattice.  **/
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > > computeSurfaceScalar(
    TriangleBoundary3D<T> &boundary, VoxelizedDomain3D<T> &voxelizedDomain,
    MultiScalarField3D<T> &scalarField, const std::vector<int> &usableFlowTypes,
    plint whichScalar = -1, bool dynamicMesh = false);

/// Requirement: particles must be of type point-particle.
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > >
    getParticlePosAndVelocity(
        MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &originalParticles);

template <typename T, template <typename U> class Descriptor>
void injectVisualScalarFieldParticles(
    MultiScalarField3D<T> &scalarField,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particleField, Box3D domain,
    Particle3D<T, Descriptor> *particleTemplate, T mostUnlikely, T mostLikely, T probability = 1.,
    plint numShotsPerCell = 1);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > >
    injectVisualScalarFieldParticles(
        MultiScalarField3D<T> &scalarField, Box3D domain,
        Particle3D<T, Descriptor> *particleTemplate, T mostUnlikely, T mostLikely,
        T probability = 1., plint numShotsPerCell = 1);

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
void scalarFieldToParticles(
    MultiScalarField3D<T> &scalar, MultiParticleField3D<ParticleFieldT> &particleField,
    Box3D domain);

template <typename T, int nDim, template <typename U> class Descriptor, class ParticleFieldT>
void tensorFieldToParticles(
    MultiTensorField3D<T, nDim> &tensor, MultiParticleField3D<ParticleFieldT> &particleField,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > > scalarFieldToSurface(
    TriangularSurfaceMesh<T> &mesh, MultiScalarField3D<T> &scalar, Box3D domain);

}  // namespace plb

#endif  // VISUAL_PARTICLE_WRAPPER_3D_H
