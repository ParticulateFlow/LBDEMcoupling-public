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

#ifndef PARTICLE_PROCESSING_FUNCTIONAL_3D_H
#define PARTICLE_PROCESSING_FUNCTIONAL_3D_H

#include <map>

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "core/functions.h"
#include "core/globalDefs.h"
#include "offLattice/triangleBoundary3D.h"

namespace plb {

/// Count the number of particles, no matter which kind, found inside the domain.
template <typename T, template <typename U> class Descriptor>
class CountParticlesFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    CountParticlesFunctional3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CountParticlesFunctional3D<T, Descriptor> *clone() const;
    plint getNumParticles() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint numParticlesId;
};

/// Count the number of particles, no matter which kind, found inside the domain.
template <typename T, template <typename U> class Descriptor>
class CountParticlesSelectiveFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    CountParticlesSelectiveFunctional3D(util::SelectInt *tags_);
    ~CountParticlesSelectiveFunctional3D();
    CountParticlesSelectiveFunctional3D(
        CountParticlesSelectiveFunctional3D<T, Descriptor> const &rhs);
    CountParticlesSelectiveFunctional3D<T, Descriptor> &operator=(
        CountParticlesSelectiveFunctional3D<T, Descriptor> const &rhs);
    void swap(CountParticlesSelectiveFunctional3D<T, Descriptor> &rhs);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CountParticlesSelectiveFunctional3D<T, Descriptor> *clone() const;
    plint getNumParticles() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint numParticlesId;
    util::SelectInt *tags;
};

/// Compute the average over all particle velocities.
template <typename T, template <typename U> class Descriptor>
class AverageParticleVelocityFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    AverageParticleVelocityFunctional3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual AverageParticleVelocityFunctional3D<T, Descriptor> *clone() const;
    Array<T, 3> getAverageParticleVelocity() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Array<plint, 3> averageVelocityId;
};

/// Copy particles of a certain tag from one field to another.
template <typename T, template <typename U> class Descriptor>
class CopySelectParticles3D : public BoxProcessingFunctional3D {
public:
    CopySelectParticles3D(util::SelectInt *tags_);
    ~CopySelectParticles3D();
    CopySelectParticles3D(CopySelectParticles3D<T, Descriptor> const &rhs);
    CopySelectParticles3D<T, Descriptor> &operator=(
        CopySelectParticles3D<T, Descriptor> const &rhs);
    void swap(CopySelectParticles3D<T, Descriptor> &rhs);
    /// Arguments: [0] From Particle-field, [1] To Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CopySelectParticles3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    util::SelectInt *tags;
};

/// Inject particles into the domain. The particles must be defined in a non-
///   parallel way, and duplicated over all processors.
template <typename T, template <typename U> class Descriptor>
class InjectParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    /// The particles are not consumed in this class. A clone of the particles is
    ///   automatically made as they are added into the domain.
    InjectParticlesFunctional3D(std::vector<Particle3D<T, Descriptor> *> &particles_);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InjectParticlesFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    std::vector<Particle3D<T, Descriptor> *> &particles;
};

/// Generate a random number of particles inside the domain. Each cell generates
///   at most one particle, with a given probability and at a random position inside
///   the cell. All particles are identical clones (except for their position).
template <typename T, template <typename U> class Descriptor>
class InjectRandomParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    InjectRandomParticlesFunctional3D(
        Particle3D<T, Descriptor> *particleTemplate_, T probabilityPerCell_);
    InjectRandomParticlesFunctional3D(InjectRandomParticlesFunctional3D<T, Descriptor> const &rhs);
    InjectRandomParticlesFunctional3D<T, Descriptor> &operator=(
        InjectRandomParticlesFunctional3D<T, Descriptor> const &rhs);
    void swap(InjectRandomParticlesFunctional3D<T, Descriptor> &rhs);
    ~InjectRandomParticlesFunctional3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InjectRandomParticlesFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
};

/// Same as the InjectRandomParticlesFunctional3D but using a parallel pseudo random
///   number generator (PPRNG) for correctness and reproducibility.
/// The rationale for the seed having a pointer type is the following. The only way to
///   control the random number generator is through the seed value. Imagine that this
///   data processor needs to be integrated and executed at every iteration with a
///   different seed each time, the value of which is controlled by the caller. A nice
///   way to achieve this is by passing the address of a variable (which is "external" to
///   the data processor) so the caller can change its value as they wish. This is why the
///   "seed" member of this class has a pointer type.
template <typename T, template <typename U> class Descriptor>
class InjectRandomParticlesFunctionalPPRNG3D : public BoxProcessingFunctional3D {
public:
    InjectRandomParticlesFunctionalPPRNG3D(
        Particle3D<T, Descriptor> *particleTemplate_, T probabilityPerCell_,
        Box3D const &boundingBox_, uint32_t const *seed_);
    InjectRandomParticlesFunctionalPPRNG3D(
        InjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> const &rhs);
    InjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> &operator=(
        InjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> const &rhs);
    void swap(InjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> &rhs);
    ~InjectRandomParticlesFunctionalPPRNG3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    plint nY, nZ;
    uint32_t const *seed;
};

/// Generate a random number of particles inside the domain. Each cell generates
///   at most one particle, with a given probability and at a random position inside
///   the cell. All particles are identical clones (except for their position).
///   Additionally to InjectRandomParticlesFunctional3D, this functional uses a bit-
///   mask to decide where to inject.
template <typename T, template <typename U> class Descriptor>
class MaskedInjectRandomParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    MaskedInjectRandomParticlesFunctional3D(
        Particle3D<T, Descriptor> *particleTemplate_, T probabilityPerCell_, int flag_);
    MaskedInjectRandomParticlesFunctional3D(
        MaskedInjectRandomParticlesFunctional3D<T, Descriptor> const &rhs);
    MaskedInjectRandomParticlesFunctional3D<T, Descriptor> &operator=(
        MaskedInjectRandomParticlesFunctional3D<T, Descriptor> const &rhs);
    void swap(MaskedInjectRandomParticlesFunctional3D<T, Descriptor> &rhs);
    ~MaskedInjectRandomParticlesFunctional3D();
    /// Arguments: [0] Particle-field; [1] Mask (int scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual MaskedInjectRandomParticlesFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    int flag;
};

/// Same as the MaskedInjectRandomParticlesFunctional3D but using a parallel pseudo random
///   number generator (PPRNG) for correctness and reproducibility.
/// The rationale for the seed having a pointer type is the following. The only way to
///   control the random number generator is through the seed value. Imagine that this
///   data processor needs to be integrated and executed at every iteration with a
///   different seed each time, the value of which is controlled by the caller. A nice
///   way to achieve this is by passing the address of a variable (which is "external" to
///   the data processor) so the caller can change its value as they wish. This is why the
///   "seed" member of this class has a pointer type.
template <typename T, template <typename U> class Descriptor>
class MaskedInjectRandomParticlesFunctionalPPRNG3D : public BoxProcessingFunctional3D {
public:
    MaskedInjectRandomParticlesFunctionalPPRNG3D(
        Particle3D<T, Descriptor> *particleTemplate_, T probabilityPerCell_, int flag_,
        Box3D const &boundingBox_, uint32_t const *seed_);
    MaskedInjectRandomParticlesFunctionalPPRNG3D(
        MaskedInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> const &rhs);
    MaskedInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> &operator=(
        MaskedInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> const &rhs);
    void swap(MaskedInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> &rhs);
    ~MaskedInjectRandomParticlesFunctionalPPRNG3D();
    /// Arguments: [0] Particle-field; [1] Mask (int scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual MaskedInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    int flag;
    plint nY, nZ;
    uint32_t const *seed;
};

template <typename T, template <typename U> class Descriptor>
class N_MaskedInjectRandomParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    N_MaskedInjectRandomParticlesFunctional3D(
        Particle3D<T, Descriptor> *particleTemplate_, T probabilityPerCell_, int flag_);
    N_MaskedInjectRandomParticlesFunctional3D(
        N_MaskedInjectRandomParticlesFunctional3D<T, Descriptor> const &rhs);
    N_MaskedInjectRandomParticlesFunctional3D<T, Descriptor> &operator=(
        N_MaskedInjectRandomParticlesFunctional3D<T, Descriptor> const &rhs);
    void swap(N_MaskedInjectRandomParticlesFunctional3D<T, Descriptor> &rhs);
    ~N_MaskedInjectRandomParticlesFunctional3D();
    /// Arguments: [0] Particle-field; [1] Mask (1D int n-tensor-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual N_MaskedInjectRandomParticlesFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    int flag;
};

/// Same as the N_MaskedInjectRandomParticlesFunctional3D but using a parallel pseudo random
///   number generator (PPRNG) for correctness and reproducibility.
/// The rationale for the seed having a pointer type is the following. The only way to
///   control the random number generator is through the seed value. Imagine that this
///   data processor needs to be integrated and executed at every iteration with a
///   different seed each time, the value of which is controlled by the caller. A nice
///   way to achieve this is by passing the address of a variable (which is "external" to
///   the data processor) so the caller can change its value as they wish. This is why the
///   "seed" member of this class has a pointer type.
template <typename T, template <typename U> class Descriptor>
class N_MaskedInjectRandomParticlesFunctionalPPRNG3D : public BoxProcessingFunctional3D {
public:
    N_MaskedInjectRandomParticlesFunctionalPPRNG3D(
        Particle3D<T, Descriptor> *particleTemplate_, T probabilityPerCell_, int flag_,
        Box3D const &boundingBox_, uint32_t const *seed_);
    N_MaskedInjectRandomParticlesFunctionalPPRNG3D(
        N_MaskedInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> const &rhs);
    N_MaskedInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> &operator=(
        N_MaskedInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> const &rhs);
    void swap(N_MaskedInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> &rhs);
    ~N_MaskedInjectRandomParticlesFunctionalPPRNG3D();
    /// Arguments: [0] Particle-field; [1] Mask (1D int n-tensor-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual N_MaskedInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    int flag;
    plint nY, nZ;
    uint32_t const *seed;
};

/// Generate a random number of point-particles inside the domain. Each cell generates
///   at most one particle, with a given probability and at a random position inside
///   the cell.
template <typename T, template <typename U> class Descriptor, class DomainFunctional>
class AnalyticalInjectRandomParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    AnalyticalInjectRandomParticlesFunctional3D(
        Particle3D<T, Descriptor> *particleTemplate_, T probabilityPerCell_,
        DomainFunctional functional_);
    AnalyticalInjectRandomParticlesFunctional3D(
        AnalyticalInjectRandomParticlesFunctional3D<T, Descriptor, DomainFunctional> const &rhs);
    AnalyticalInjectRandomParticlesFunctional3D<T, Descriptor, DomainFunctional> &operator=(
        AnalyticalInjectRandomParticlesFunctional3D<T, Descriptor, DomainFunctional> const &rhs);
    void swap(AnalyticalInjectRandomParticlesFunctional3D<T, Descriptor, DomainFunctional> &rhs);
    ~AnalyticalInjectRandomParticlesFunctional3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual AnalyticalInjectRandomParticlesFunctional3D<T, Descriptor, DomainFunctional> *clone()
        const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    DomainFunctional functional;
};

/// Same as the AnalyticalInjectRandomParticlesFunctional3D but using a parallel pseudo random
///   number generator (PPRNG) for correctness and reproducibility.
/// The rationale for the seed having a pointer type is the following. The only way to
///   control the random number generator is through the seed value. Imagine that this
///   data processor needs to be integrated and executed at every iteration with a
///   different seed each time, the value of which is controlled by the caller. A nice
///   way to achieve this is by passing the address of a variable (which is "external" to
///   the data processor) so the caller can change its value as they wish. This is why the
///   "seed" member of this class has a pointer type.
template <typename T, template <typename U> class Descriptor, class DomainFunctional>
class AnalyticalInjectRandomParticlesFunctionalPPRNG3D : public BoxProcessingFunctional3D {
public:
    AnalyticalInjectRandomParticlesFunctionalPPRNG3D(
        Particle3D<T, Descriptor> *particleTemplate_, T probabilityPerCell_,
        DomainFunctional functional_, Box3D const &boundingBox_, uint32_t const *seed_);
    AnalyticalInjectRandomParticlesFunctionalPPRNG3D(
        AnalyticalInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor, DomainFunctional> const
            &rhs);
    AnalyticalInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor, DomainFunctional> &operator=(
        AnalyticalInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor, DomainFunctional> const
            &rhs);
    void swap(
        AnalyticalInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor, DomainFunctional> &rhs);
    ~AnalyticalInjectRandomParticlesFunctionalPPRNG3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual AnalyticalInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor, DomainFunctional>
        *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    DomainFunctional functional;
    plint nY, nZ;
    uint32_t const *seed;
};

/// Generate a random number of point-particles inside the domain. Each cell generates
///   at most one particle, with a given probability and at a random position inside
///   the cell. Additionally to analytically-inject, this functional uses a bit-
///   mask to decide where to inject.
template <typename T, template <typename U> class Descriptor, class DomainFunctional>
class MaskedAnalyticalInjectRandomParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    MaskedAnalyticalInjectRandomParticlesFunctional3D(
        Particle3D<T, Descriptor> *particleTemplate_, T probabilityPerCell_,
        DomainFunctional functional_, int flag_);
    MaskedAnalyticalInjectRandomParticlesFunctional3D(
        MaskedAnalyticalInjectRandomParticlesFunctional3D<T, Descriptor, DomainFunctional> const
            &rhs);
    MaskedAnalyticalInjectRandomParticlesFunctional3D<T, Descriptor, DomainFunctional> &operator=(
        MaskedAnalyticalInjectRandomParticlesFunctional3D<T, Descriptor, DomainFunctional> const
            &rhs);
    void swap(
        MaskedAnalyticalInjectRandomParticlesFunctional3D<T, Descriptor, DomainFunctional> &rhs);
    ~MaskedAnalyticalInjectRandomParticlesFunctional3D();
    /// Arguments: Particle-field, Mask.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual MaskedAnalyticalInjectRandomParticlesFunctional3D<T, Descriptor, DomainFunctional>
        *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    DomainFunctional functional;
    int flag;
};

/// Same as the MaskedAnalyticalInjectRandomParticlesFunctional3D but using a parallel pseudo random
///   number generator (PPRNG) for correctness and reproducibility.
/// The rationale for the seed having a pointer type is the following. The only way to
///   control the random number generator is through the seed value. Imagine that this
///   data processor needs to be integrated and executed at every iteration with a
///   different seed each time, the value of which is controlled by the caller. A nice
///   way to achieve this is by passing the address of a variable (which is "external" to
///   the data processor) so the caller can change its value as they wish. This is why the
///   "seed" member of this class has a pointer type.
template <typename T, template <typename U> class Descriptor, class DomainFunctional>
class MaskedAnalyticalInjectRandomParticlesFunctionalPPRNG3D : public BoxProcessingFunctional3D {
public:
    MaskedAnalyticalInjectRandomParticlesFunctionalPPRNG3D(
        Particle3D<T, Descriptor> *particleTemplate_, T probabilityPerCell_,
        DomainFunctional functional_, int flag_, Box3D const &boundingBox_, uint32_t const *seed_);
    MaskedAnalyticalInjectRandomParticlesFunctionalPPRNG3D(
        MaskedAnalyticalInjectRandomParticlesFunctionalPPRNG3D<
            T, Descriptor, DomainFunctional> const &rhs);
    MaskedAnalyticalInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor, DomainFunctional>
        &operator=(MaskedAnalyticalInjectRandomParticlesFunctionalPPRNG3D<
                   T, Descriptor, DomainFunctional> const &rhs);
    void swap(
        MaskedAnalyticalInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor, DomainFunctional>
            &rhs);
    ~MaskedAnalyticalInjectRandomParticlesFunctionalPPRNG3D();
    /// Arguments: Particle-field, Mask.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual MaskedAnalyticalInjectRandomParticlesFunctionalPPRNG3D<T, Descriptor, DomainFunctional>
        *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    DomainFunctional functional;
    int flag;
    plint nY, nZ;
    uint32_t const *seed;
};

/// Generate equally spaced particles inside each cell. Every cell generates
///   nx particles in its x-direction, ny in its y-direction and nz in its
///   z-direction (nx * ny * nz in total in each cell).
///   All particles are identical clones (except for their position).
template <typename T, template <typename U> class Descriptor>
class InjectEquallySpacedParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    InjectEquallySpacedParticlesFunctional3D(
        Particle3D<T, Descriptor> *particleTemplate_, plint nx_, plint ny_, plint nz_);
    InjectEquallySpacedParticlesFunctional3D(
        InjectEquallySpacedParticlesFunctional3D<T, Descriptor> const &rhs);
    InjectEquallySpacedParticlesFunctional3D<T, Descriptor> &operator=(
        InjectEquallySpacedParticlesFunctional3D<T, Descriptor> const &rhs);
    void swap(InjectEquallySpacedParticlesFunctional3D<T, Descriptor> &rhs);
    ~InjectEquallySpacedParticlesFunctional3D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InjectEquallySpacedParticlesFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    plint nx, ny, nz;
};

/// Generate equally spaced particles inside each cell. Every cell generates
///   nx particles in its x-direction, ny in its y-direction and nz in its
///   z-direction (nx * ny * nz in total in each cell). This functional uses
///   a bit-mask to decide where to inject. All particles are identical clones
///   (except for their position).
template <typename T, template <typename U> class Descriptor>
class MaskedInjectEquallySpacedParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    MaskedInjectEquallySpacedParticlesFunctional3D(
        Particle3D<T, Descriptor> *particleTemplate_, plint nx_, plint ny_, plint nz_, int flag_);
    MaskedInjectEquallySpacedParticlesFunctional3D(
        MaskedInjectEquallySpacedParticlesFunctional3D<T, Descriptor> const &rhs);
    MaskedInjectEquallySpacedParticlesFunctional3D<T, Descriptor> &operator=(
        MaskedInjectEquallySpacedParticlesFunctional3D<T, Descriptor> const &rhs);
    void swap(MaskedInjectEquallySpacedParticlesFunctional3D<T, Descriptor> &rhs);
    ~MaskedInjectEquallySpacedParticlesFunctional3D();
    /// Arguments: Particle-field, Mask.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual MaskedInjectEquallySpacedParticlesFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    plint nx, ny, nz;
    int flag;
};

/// Generate equally spaced particles inside each cell. Every cell generates
///   nx particles in its x-direction, ny in its y-direction and nz in its
///   z-direction (nx * ny * nz in total in each cell). This functional uses
///   a bit-mask (in the form of a MultiNTensorField3D) to decide where to inject.
///   All particles are identical clones (except for their position).
template <typename T, template <typename U> class Descriptor>
class N_MaskedInjectEquallySpacedParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    N_MaskedInjectEquallySpacedParticlesFunctional3D(
        Particle3D<T, Descriptor> *particleTemplate_, plint nx_, plint ny_, plint nz_, int flag_);
    N_MaskedInjectEquallySpacedParticlesFunctional3D(
        N_MaskedInjectEquallySpacedParticlesFunctional3D<T, Descriptor> const &rhs);
    N_MaskedInjectEquallySpacedParticlesFunctional3D<T, Descriptor> &operator=(
        N_MaskedInjectEquallySpacedParticlesFunctional3D<T, Descriptor> const &rhs);
    void swap(N_MaskedInjectEquallySpacedParticlesFunctional3D<T, Descriptor> &rhs);
    ~N_MaskedInjectEquallySpacedParticlesFunctional3D();
    /// Arguments: Particle-field, Mask.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual N_MaskedInjectEquallySpacedParticlesFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle3D<T, Descriptor> *particleTemplate;
    plint nx, ny, nz;
    int flag;
};

/// Remove all particles from a given domain.
template <typename T, template <typename U> class Descriptor>
class AbsorbParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    /// Argument: Particle-field.
    virtual AbsorbParticlesFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/// Remove all particles based on a mask.
template <typename T, template <typename U> class Descriptor>
class MaskedAbsorbParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    MaskedAbsorbParticlesFunctional3D(int whichFlag_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    /// Arguments: [0] Particle-field; [1] Mask (int scalar-field)
    virtual MaskedAbsorbParticlesFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int whichFlag;
};

/// Remove all particles from a given domain.
template <typename T, template <typename U> class Descriptor>
class AbsorbParticlesFunctionalSelective3D : public BoxProcessingFunctional3D {
public:
    AbsorbParticlesFunctionalSelective3D(util::SelectInt *tags_);
    ~AbsorbParticlesFunctionalSelective3D();
    AbsorbParticlesFunctionalSelective3D(
        AbsorbParticlesFunctionalSelective3D<T, Descriptor> const &rhs);
    AbsorbParticlesFunctionalSelective3D<T, Descriptor> &operator=(
        AbsorbParticlesFunctionalSelective3D<T, Descriptor> const &rhs);
    void swap(AbsorbParticlesFunctionalSelective3D<T, Descriptor> &rhs);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    /// Argument: Particle-field.
    virtual AbsorbParticlesFunctionalSelective3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    util::SelectInt *tags;
};

/// Find particles injected inside wall nodes and remove them.
template <typename T, template <typename U> class Descriptor>
class RemoveParticlesFromWall3D : public BoxProcessingFunctional3D {
public:
    RemoveParticlesFromWall3D(int wallFlag_);
    /// Arguments: [0] Particle-field [1] Flag-matrix
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual RemoveParticlesFromWall3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int wallFlag;  // Value that represents the wall in the flag matrix.
};

/// Find particles close to a wall and change their positions so they are pushed back to the flow
/// field.
template <typename T, template <typename U> class Descriptor>
class PushParticlesAwayFromWall3D : public BoxProcessingFunctional3D {
public:
    PushParticlesAwayFromWall3D(T cutOffValue_, T movingDistance_, int wallFlag_, int fluidFlag_);
    /// Arguments: [0] Particle-field [1] Flag-matrix
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual PushParticlesAwayFromWall3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T cutOffValue;     // When the speed of the particle drops below sqrt(cutOffValue), then this
                       // particle is a candidate for pushing.
    T movingDistance;  // This is the distance the particles will be moved.
    int wallFlag;      // Value that represents the wall nodes in the flag matrix.
    int fluidFlag;     // Value that represents the fluid nodes in the flag matrix.
};

/// Execute the particle-fluid interaction step (during which the particles
///   don't move and the fluid doesn't change).
template <typename T, template <typename U> class Descriptor>
class FluidToParticleCoupling3D : public BoxProcessingFunctional3D {
public:
    /// Particle speed = scaling*fluid speed.
    FluidToParticleCoupling3D(T scaling_);
    /// Arguments: [0] Particle-field; [1] Fluid.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual FluidToParticleCoupling3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T scaling;
};

template <typename T, template <typename U> class Descriptor>
class VelocityToParticleCoupling3D : public BoxProcessingFunctional3D {
public:
    /// Particle speed = scaling*fluid speed.
    VelocityToParticleCoupling3D(T scaling_);
    /// Arguments: [0] Particle-field; [1] Velocity.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual VelocityToParticleCoupling3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T scaling;
};

template <typename T, template <typename U> class Descriptor>
class N_VelocityToParticleCoupling3D : public BoxProcessingFunctional3D {
public:
    /// Particle speed = scaling*fluid speed.
    N_VelocityToParticleCoupling3D(T scaling_);
    /// Arguments: [0] Particle-field; [1] Velocity (scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual N_VelocityToParticleCoupling3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T scaling;
};

template <typename T, template <typename U> class Descriptor>
class RhoBarJtoParticleCoupling3D : public BoxProcessingFunctional3D {
public:
    /// Particle speed = scaling*fluid speed.
    RhoBarJtoParticleCoupling3D(bool velIsJ_, T scaling_);
    /// Arguments: [0] Particle-field; [1] rhoBarJ.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual RhoBarJtoParticleCoupling3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    bool velIsJ;
    T scaling;
};

/// Execute the iteration step during which particles advance.
/// WARNING: This data processor does not update the envelope of the multi-particle
/// field after its execution.
template <typename T, template <typename U> class Descriptor>
class AdvanceParticlesFunctional3D : public BoxProcessingFunctional3D {
public:
    /// When the speed of a particle drops below sqrt(cutOffValue),
    ///   the particle is eliminated. Negative cutOffValue means no cutoff.
    AdvanceParticlesFunctional3D(T cutOffValue_ = -1.);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual AdvanceParticlesFunctional3D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T cutOffValue;
};

/// Execute the iteration step during which particles advance, on the whole domain
/** The data processor's domain indication is being ignored. This works also with periodicity. **/
template <typename T, template <typename U> class Descriptor>
class AdvanceParticlesEveryWhereFunctional3D : public BoxProcessingFunctional3D {
public:
    /// When the speed of a particle drops below sqrt(cutOffValue),
    ///   the particle is eliminated. Negative cutOffValue means no cutoff.
    AdvanceParticlesEveryWhereFunctional3D(T cutOffValue_ = -1.);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual AdvanceParticlesEveryWhereFunctional3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T cutOffValue;
};

/* ******** VerletUpdateVelocity3D *********************************** */

/// Update the velocity to complete an iteration of the Verlet algorithm. Works
/// with Verlet particles only.
template <typename T, template <typename U> class Descriptor>
class VerletUpdateVelocity3D : public BoxProcessingFunctional3D {
public:
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual VerletUpdateVelocity3D<T, Descriptor> *clone() const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    bool projectForce;
    Array<T, 3> planeNormal;
};

/* ******** VerletUpdateVelocitySelective3D *********************************** */

/// Update the velocity to complete an iteration of the Verlet algorithm. Works
/// with Verlet particles only. Acts only on particles with the specified tag.
template <typename T, template <typename U> class Descriptor>
class VerletUpdateVelocitySelective3D : public BoxProcessingFunctional3D {
public:
    VerletUpdateVelocitySelective3D(util::SelectInt *tags_);
    ~VerletUpdateVelocitySelective3D();
    VerletUpdateVelocitySelective3D(VerletUpdateVelocitySelective3D<T, Descriptor> const &rhs);
    VerletUpdateVelocitySelective3D<T, Descriptor> &operator=(
        VerletUpdateVelocitySelective3D<T, Descriptor> const &rhs);
    void swap(VerletUpdateVelocitySelective3D<T, Descriptor> &rhs);
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual VerletUpdateVelocitySelective3D<T, Descriptor> *clone() const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    util::SelectInt *tags;
};

/* ******** VerletParticleInteractionForce3D *********************************** */

// This data processor defines the pair-wise interaction force between Verlet particles.
// The cutOffLength is the cut-off length for particle interaction. One needs to make sure
// the interaction force is negligible at this distance. The halfWidth is the width of
// communication envelope needed to cope with neighborhoods within the cut-off length.
// The particle-field used in this data processor, must have an envelope width at least
// as large as the halfWidth provided. The algorithm is not applied to particles that have
// a tag equal to the exclusionTag (if provided).
template <typename T, template <typename U> class Descriptor>
class VerletParticleInteractionForce3D : public BoxProcessingFunctional3D {
public:
    VerletParticleInteractionForce3D(
        T forceAmplitude_, Array<T, 3> const &bodyAcceleration_, T cutOffLength_, plint halfWidth_,
        plint exclusionTag_);
    VerletParticleInteractionForce3D(
        T forceAmplitude_, Array<T, 3> const &bodyAcceleration_, T cutOffLength_, plint halfWidth_);
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual VerletParticleInteractionForce3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T forceAmplitude;
    Array<T, 3> bodyAcceleration;
    T cutOffLength;
    plint halfWidth;
    plint exclusionTag;
    bool hasExclusionTag;
};

template <typename T, template <typename U> class Descriptor>
void addWallParticles(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles,
    TriangleBoundary3D<T> &boundary);

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
void addWallParticlesGeneric(
    MultiParticleField3D<ParticleFieldT> &particles, TriangleBoundary3D<T> &boundary);

/// Count the number of particles at each cell node and add the result to the scalar field.
template <typename T, template <typename U> class Descriptor>
class CountAndAccumulateParticles3D : public BoxProcessingFunctional3D {
public:
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CountAndAccumulateParticles3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/// Count the number of particles with a given tag at each cell node and add the result to the
/// scalar field.
template <typename T, template <typename U> class Descriptor>
class CountAndAccumulateTaggedParticles3D : public BoxProcessingFunctional3D {
public:
    CountAndAccumulateTaggedParticles3D(plint tag_);
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CountAndAccumulateTaggedParticles3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint tag;
};

/// Count the number of particles (with a given tag) at each refined cell node, and add the result
/// to the scalar field which is refined (defined on a refined grid with respect to the particle
/// grid). The particles which belong to each "sub-volume" of the refined scalar grid contained in
/// the "big volume" of the particle grid, must be identified, counted and accumulated.
template <typename T, template <typename U> class Descriptor>
class CountAndAccumulateTaggedParticlesRefined3D : public BoxProcessingFunctional3D {
public:
    CountAndAccumulateTaggedParticlesRefined3D(plint tag_, plint dxScale_);
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box3D coarseDomain, std::vector<AtomicBlock3D *> fields);
    virtual CountAndAccumulateTaggedParticlesRefined3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    /// These helper functions are also implemented in the ParticleField3D class,
    /// but we need to re-implement them here, since we need them for the refined
    /// scalar field.
    static plint nearestCell(T pos);
    static void computeGridPosition(
        Array<T, 3> const &position, Dot3D const &location, plint &iX, plint &iY, plint &iZ);

private:
    plint tag;
    plint dxScale;
};

/// Count the number of particles with given tags at each cell node and place the result to the
/// scalar field.
template <typename T, template <typename U> class Descriptor>
class CountTaggedParticles3D : public BoxProcessingFunctional3D {
public:
    CountTaggedParticles3D(util::SelectInt *tags_);
    ~CountTaggedParticles3D();
    CountTaggedParticles3D(CountTaggedParticles3D<T, Descriptor> const &rhs);
    CountTaggedParticles3D<T, Descriptor> &operator=(
        CountTaggedParticles3D<T, Descriptor> const &rhs);
    void swap(CountTaggedParticles3D<T, Descriptor> &rhs);
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual CountTaggedParticles3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    util::SelectInt *tags;
};

/// Count the number of particles with given tags inside a neighborhood of each cell node,
/// compute the corresponding concentration value and place the result in the scalar field.
/// The parameter halfWidth must be at most equal to the envelope width of the particle field.
template <typename T, template <typename U> class Descriptor>
class ComputeConcentrationOfTaggedParticles3D : public BoxProcessingFunctional3D {
public:
    ComputeConcentrationOfTaggedParticles3D(util::SelectInt *tags_, plint halfWidth_);
    ~ComputeConcentrationOfTaggedParticles3D();
    ComputeConcentrationOfTaggedParticles3D(
        ComputeConcentrationOfTaggedParticles3D<T, Descriptor> const &rhs);
    ComputeConcentrationOfTaggedParticles3D<T, Descriptor> &operator=(
        ComputeConcentrationOfTaggedParticles3D<T, Descriptor> const &rhs);
    void swap(ComputeConcentrationOfTaggedParticles3D<T, Descriptor> &rhs);
    /// Arguments: [0] Particle-field; [1] Concentration (real-valued scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ComputeConcentrationOfTaggedParticles3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    util::SelectInt *tags;
    plint halfWidth;
};

/// Count the number of particles with given tags inside a neighborhood of each cell node,
/// compute the corresponding concentration value and place the result in the scalar field,
/// only for cells that have a specific mask value.
/// The parameter halfWidth must be at most equal to the envelope width of the particle field.
template <typename T, template <typename U> class Descriptor>
class MaskedComputeConcentrationOfTaggedParticles3D : public BoxProcessingFunctional3D {
public:
    MaskedComputeConcentrationOfTaggedParticles3D(
        util::SelectInt *tags_, plint halfWidth_, int whichFlag_);
    ~MaskedComputeConcentrationOfTaggedParticles3D();
    MaskedComputeConcentrationOfTaggedParticles3D(
        MaskedComputeConcentrationOfTaggedParticles3D<T, Descriptor> const &rhs);
    MaskedComputeConcentrationOfTaggedParticles3D<T, Descriptor> &operator=(
        MaskedComputeConcentrationOfTaggedParticles3D<T, Descriptor> const &rhs);
    void swap(MaskedComputeConcentrationOfTaggedParticles3D<T, Descriptor> &rhs);
    /// Arguments: [0] Particle-field; [1] Mask (int scalar-field); [2] Concentration (real-valued
    /// scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual MaskedComputeConcentrationOfTaggedParticles3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    util::SelectInt *tags;
    plint halfWidth;
    int whichFlag;
};

template <typename T, template <typename U> class Descriptor>
class N_MaskedComputeConcentrationOfTaggedParticles3D : public BoxProcessingFunctional3D {
public:
    N_MaskedComputeConcentrationOfTaggedParticles3D(
        util::SelectInt *tags_, plint halfWidth_, int whichFlag_);
    ~N_MaskedComputeConcentrationOfTaggedParticles3D();
    N_MaskedComputeConcentrationOfTaggedParticles3D(
        N_MaskedComputeConcentrationOfTaggedParticles3D<T, Descriptor> const &rhs);
    N_MaskedComputeConcentrationOfTaggedParticles3D<T, Descriptor> &operator=(
        N_MaskedComputeConcentrationOfTaggedParticles3D<T, Descriptor> const &rhs);
    void swap(N_MaskedComputeConcentrationOfTaggedParticles3D<T, Descriptor> &rhs);
    /// Arguments: [0] Particle-field; [1] Mask (1D int n-tensor-field); [2] Concentration
    /// (real-valued scalar-field).
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual N_MaskedComputeConcentrationOfTaggedParticles3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    util::SelectInt *tags;
    plint halfWidth;
    int whichFlag;
};

template <
    typename T, template <typename U> class Descriptor,
    template <typename T_, template <typename U_> class Descriptor_> class ParticleFieldT>
plint countParticles(
    MultiParticleField3D<ParticleFieldT<T, Descriptor> > &particles, Box3D const &domain);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T_, template <typename U_> class Descriptor_> class ParticleFieldT>
plint countParticles(
    MultiParticleField3D<ParticleFieldT<T, Descriptor> > &particles, Box3D const &domain,
    util::SelectInt *tags);

template <typename T, template <typename U> class Descriptor>
void injectParticles(
    std::vector<Particle3D<T, Descriptor> *> &injectedParticles,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, Box3D domain);

template <typename T, template <typename U> class Descriptor, class ParticleFieldT>
void injectParticles(
    std::vector<Particle3D<T, Descriptor> *> &injectedParticles,
    MultiParticleField3D<ParticleFieldT> &particles, Box3D domain);

/* Iterations of a passive-scalar fluid-particle system:
 * =====================================================
 *
 * Note: The difficulty comes from the fact that particle-fields may have a larger
 *   envelope than the fluid. When advancing particles on bulk and envelope, the
 *   velocity data from the fluid is therefore not necessarily locally available.
 *   The velocity is therefore first stored in the particle (in the bulk), and then
 *   communicated to the envelopes.
 *
 * --- Particles are at time t, fluid is at time t, defined on bulk and envelope. ---
 *  1. Fluid collideAndStream().
 *  2. Particle advance (bulk+envelope). ==> Particles at time t on bulk (needs no communication).
 *  3. Fluid communication ==> Fluid at time t+1.
 *  4. Particle interact (bulk domain) with velocity at time t+1.
 *  5. Particle communication ==> Particle at time t+1.
 */

}  // namespace plb

#endif  // PARTICLE_PROCESSING_FUNCTIONAL_3D_H
