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

#ifndef PARTICLE_PROCESSING_FUNCTIONAL_2D_H
#define PARTICLE_PROCESSING_FUNCTIONAL_2D_H

#include <map>

#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "core/globalDefs.h"

namespace plb {

/// Count the number of particles, no matter which kind, found inside the domain.
template <typename T, template <typename U> class Descriptor>
class CountParticlesFunctional2D : public PlainReductiveBoxProcessingFunctional2D {
public:
    CountParticlesFunctional2D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual CountParticlesFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    plint getNumParticles() const;

private:
    plint numParticlesId;
};

/// Count the number of particles, no matter which kind, found inside the domain.
template <typename T, template <typename U> class Descriptor>
class CountParticlesSelectiveFunctional2D : public PlainReductiveBoxProcessingFunctional2D {
public:
    CountParticlesSelectiveFunctional2D(util::SelectInt *tags_);
    ~CountParticlesSelectiveFunctional2D();
    CountParticlesSelectiveFunctional2D(
        CountParticlesSelectiveFunctional2D<T, Descriptor> const &rhs);
    CountParticlesSelectiveFunctional2D<T, Descriptor> &operator=(
        CountParticlesSelectiveFunctional2D<T, Descriptor> const &rhs);
    void swap(CountParticlesSelectiveFunctional2D<T, Descriptor> &rhs);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual CountParticlesSelectiveFunctional2D<T, Descriptor> *clone() const;
    plint getNumParticles() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint numParticlesId;
    util::SelectInt *tags;
};

/// Compute the average over all particle velocities.
template <typename T, template <typename U> class Descriptor>
class AverageParticleVelocityFunctional2D : public PlainReductiveBoxProcessingFunctional2D {
public:
    AverageParticleVelocityFunctional2D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual AverageParticleVelocityFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    Array<T, 2> getAverageParticleVelocity() const;

private:
    Array<plint, 2> averageVelocityId;
};

/// Inject particles into the domain. The particles must be defined in a non-
///   parallel way, and duplicated over all processors.
template <typename T, template <typename U> class Descriptor>
class InjectParticlesFunctional2D : public BoxProcessingFunctional2D {
public:
    /// The particles are not consumed in this class. A clone of the particles is
    ///   automatically made as they are added into the domain.
    InjectParticlesFunctional2D(std::vector<Particle2D<T, Descriptor> *> &particles_);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual InjectParticlesFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    std::vector<Particle2D<T, Descriptor> *> &particles;
};

/// Generate a random number of particles inside the domain. Each cell generates
///   at most one cell, with a given probability and at a random position inside
///   the cell.
template <typename T, template <typename U> class Descriptor>
class InjectRandomPointParticlesFunctional2D : public BoxProcessingFunctional2D {
public:
    InjectRandomPointParticlesFunctional2D(plint tag_, T probabilityPerCell_);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual InjectRandomPointParticlesFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint tag;
    T probabilityPerCell;
};

/// Same as the InjectRandomPointParticlesFunctional2D but using a parallel pseudo random
///   number generator (PPRNG) for correctness and reproducibility.
/// The rationale for the seed having a pointer type is the following. The only way to
///   control the random number generator is through the seed value. Imagine that this
///   data processor needs to be integrated and executed at every iteration with a
///   different seed each time, the value of which is controlled by the caller. A nice
///   way to achieve this is by passing the address of a variable (which is "external" to
///   the data processor) so the caller can change its value as they wish. This is why the
///   "seed" member of this class has a pointer type.
template <typename T, template <typename U> class Descriptor>
class InjectRandomPointParticlesFunctionalPPRNG2D : public BoxProcessingFunctional2D {
public:
    InjectRandomPointParticlesFunctionalPPRNG2D(
        plint tag_, T probabilityPerCell_, Box2D const &boundingBox_, uint32_t const *seed_);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual InjectRandomPointParticlesFunctionalPPRNG2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint tag;
    T probabilityPerCell;
    plint nY;
    uint32_t const *seed;
};

/// Generate a random number of point-particles inside the domain. Each cell generates
///   at most one particle, with a given probability and at a random position inside
///   the cell.
template <typename T, template <typename U> class Descriptor, class DomainFunctional>
class AnalyticalInjectRandomParticlesFunctional2D : public BoxProcessingFunctional2D {
public:
    AnalyticalInjectRandomParticlesFunctional2D(
        Particle2D<T, Descriptor> *particleTemplate_, T probabilityPerCell_,
        DomainFunctional functional_);
    AnalyticalInjectRandomParticlesFunctional2D(
        AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional> const &rhs);
    AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional> &operator=(
        AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional> const &rhs);
    void swap(AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional> &rhs);
    ~AnalyticalInjectRandomParticlesFunctional2D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional> *clone()
        const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle2D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    DomainFunctional functional;
};

/// Same as the AnalyticalInjectRandomParticlesFunctional2D but using a parallel pseudo random
///   number generator (PPRNG) for correctness and reproducibility.
/// The rationale for the seed having a pointer type is the following. The only way to
///   control the random number generator is through the seed value. Imagine that this
///   data processor needs to be integrated and executed at every iteration with a
///   different seed each time, the value of which is controlled by the caller. A nice
///   way to achieve this is by passing the address of a variable (which is "external" to
///   the data processor) so the caller can change its value as they wish. This is why the
///   "seed" member of this class has a pointer type.
template <typename T, template <typename U> class Descriptor, class DomainFunctional>
class AnalyticalInjectRandomParticlesFunctionalPPRNG2D : public BoxProcessingFunctional2D {
public:
    AnalyticalInjectRandomParticlesFunctionalPPRNG2D(
        Particle2D<T, Descriptor> *particleTemplate_, T probabilityPerCell_,
        DomainFunctional functional_, Box2D const &boundingBox_, uint32_t const *seed_);
    AnalyticalInjectRandomParticlesFunctionalPPRNG2D(
        AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional> const
            &rhs);
    AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional> &operator=(
        AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional> const
            &rhs);
    void swap(
        AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional> &rhs);
    ~AnalyticalInjectRandomParticlesFunctionalPPRNG2D();
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>
        *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Particle2D<T, Descriptor> *particleTemplate;
    T probabilityPerCell;
    DomainFunctional functional;
    plint nY;
    uint32_t const *seed;
};

/// Remove all particles from a given domain.
template <typename T, template <typename U> class Descriptor>
class AbsorbParticlesFunctional2D : public BoxProcessingFunctional2D {
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    /// Argument: Particle-field.
    virtual AbsorbParticlesFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/// Execute the particle-fluid interaction step (during which the particles
///   don't move and the fluid doesn't change).
template <typename T, template <typename U> class Descriptor>
class FluidToParticleCoupling2D : public BoxProcessingFunctional2D {
public:
    /// Particle speed = scaling*fluid speed.
    FluidToParticleCoupling2D(T scaling_);
    /// Arguments: [0] Particle-field; [1] Fluid.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual FluidToParticleCoupling2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T scaling;
};

template <typename T, template <typename U> class Descriptor>
class VelocityToParticleCoupling2D : public BoxProcessingFunctional2D {
public:
    /// Particle speed = scaling*fluid speed.
    VelocityToParticleCoupling2D(T scaling_);
    /// Arguments: [0] Particle-field; [1] Velocity.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual VelocityToParticleCoupling2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T scaling;
};

template <typename T, template <typename U> class Descriptor>
class N_VelocityToParticleCoupling2D : public BoxProcessingFunctional2D {
public:
    /// Particle speed = scaling*fluid speed.
    N_VelocityToParticleCoupling2D(T scaling_);
    /// Arguments: [0] Particle-field; [1] Velocity.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual N_VelocityToParticleCoupling2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T scaling;
};

template <typename T, template <typename U> class Descriptor>
class RhoBarJtoParticleCoupling2D : public BoxProcessingFunctional2D {
public:
    /// Particle speed = scaling*fluid speed.
    RhoBarJtoParticleCoupling2D(bool velIsJ_, T scaling_);
    /// Arguments: [0] Particle-field; [1] rhoBarJ.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual RhoBarJtoParticleCoupling2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    bool velIsJ;
    T scaling;
};

/// Execute the iteration step during which particles advance.
template <typename T, template <typename U> class Descriptor>
class AdvanceParticlesFunctional2D : public BoxProcessingFunctional2D {
public:
    /// When the speed of a particle drops below sqrt(cutOffValue),
    ///   the particle is eliminated. Negative cutOffValue means no cutoff.
    AdvanceParticlesFunctional2D(T cutOffValue_ = -1.);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual AdvanceParticlesFunctional2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T cutOffValue;
};

/// Execute the iteration step during which particles advance, on the whole domain
/** The data processor's domain indication is being ignored. This works also with periodicity. **/
template <typename T, template <typename U> class Descriptor>
class AdvanceParticlesEveryWhereFunctional2D : public BoxProcessingFunctional2D {
public:
    /// When the speed of a particle drops below sqrt(cutOffValue),
    ///   the particle is eliminated. Negative cutOffValue means no cutoff.
    AdvanceParticlesEveryWhereFunctional2D(T cutOffValue_ = -1.);
    /// Argument: Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual AdvanceParticlesEveryWhereFunctional2D<T, Descriptor> *clone() const;
    //     virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T cutOffValue;
};

/* ******** VerletUpdateVelocity2D *********************************** */

/// Update the velocity to complete an iteration of the Verlet algorithm. Works
/// with Verlet particles only.
template <typename T, template <typename U> class Descriptor>
class VerletUpdateVelocity2D : public BoxProcessingFunctional2D {
public:
    /// Arguments: [0] Particle-field.
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual VerletUpdateVelocity2D<T, Descriptor> *clone() const;
    virtual void getModificationPattern(std::vector<bool> &isWritten) const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    bool projectForce;
    Array<T, 2> planeNormal;
};

/// Count the number of particles at each cell node and add the result to the scalar field.
template <typename T, template <typename U> class Descriptor>
class CountAndAccumulateParticles2D : public BoxProcessingFunctional2D {
public:
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual CountAndAccumulateParticles2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

/// Count the number of particles with a given tag at each cell node and add the result to the
/// scalar field.
template <typename T, template <typename U> class Descriptor>
class CountAndAccumulateTaggedParticles2D : public BoxProcessingFunctional2D {
public:
    CountAndAccumulateTaggedParticles2D(plint tag_);
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual CountAndAccumulateTaggedParticles2D<T, Descriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    plint tag;
};

/// Count the number of particles with given tags at each cell node and place the result to the
/// scalar field.
template <typename T, template <typename U> class Descriptor>
class CountTaggedParticles2D : public BoxProcessingFunctional2D {
public:
    CountTaggedParticles2D(util::SelectInt *tags_);
    ~CountTaggedParticles2D();
    CountTaggedParticles2D(CountTaggedParticles2D<T, Descriptor> const &rhs);
    CountTaggedParticles2D<T, Descriptor> &operator=(
        CountTaggedParticles2D<T, Descriptor> const &rhs);
    void swap(CountTaggedParticles2D<T, Descriptor> &rhs);
    /// Arguments: [0] Particle-field; [1] Number of particles (plint scalar-field).
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual CountTaggedParticles2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    util::SelectInt *tags;
};

template <
    typename T, template <typename U> class Descriptor,
    template <typename T_, template <typename U_> class Descriptor_> class ParticleFieldT>
plint countParticles(
    MultiParticleField2D<ParticleFieldT<T, Descriptor> > &particles, Box2D const &domain);

template <
    typename T, template <typename U> class Descriptor,
    template <typename T_, template <typename U_> class Descriptor_> class ParticleFieldT>
plint countParticles(
    MultiParticleField2D<ParticleFieldT<T, Descriptor> > &particles, Box2D const &domain,
    util::SelectInt *tags);

template <typename T, template <typename U> class Descriptor>
void injectParticles(
    std::vector<Particle2D<T, Descriptor> *> &injectedParticles,
    MultiParticleField2D<DenseParticleField2D<T, Descriptor> > &particles, Box2D domain);

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

#endif  // PARTICLE_PROCESSING_FUNCTIONAL_2D_H
