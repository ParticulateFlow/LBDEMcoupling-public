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

#ifndef PARTICLE_PROCESSING_FUNCTIONAL_2D_HH
#define PARTICLE_PROCESSING_FUNCTIONAL_2D_HH

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "dataProcessors/metaStuffFunctional2D.h"
#include "particles/particleField2D.h"
#include "particles/particleProcessingFunctional2D.h"
#include "sitmo/prng_engine.hpp"

namespace plb {

/* ******** CountParticlesFunctional2D *********************************** */

template <typename T, template <typename U> class Descriptor>
CountParticlesFunctional2D<T, Descriptor>::CountParticlesFunctional2D() :
    numParticlesId(this->getStatistics().subscribeIntSum())
{ }

template <typename T, template <typename U> class Descriptor>
void CountParticlesFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    std::vector<Particle2D<T, Descriptor> *> particles;
    particleField.findParticles(domain, particles);
    this->getStatistics().gatherIntSum(numParticlesId, (plint)particles.size());
}

template <typename T, template <typename U> class Descriptor>
CountParticlesFunctional2D<T, Descriptor> *CountParticlesFunctional2D<T, Descriptor>::clone() const
{
    return new CountParticlesFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CountParticlesFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT CountParticlesFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
plint CountParticlesFunctional2D<T, Descriptor>::getNumParticles() const
{
    return this->getStatistics().getIntSum(numParticlesId);
}

/* ******** CountParticlesSelectiveFunctional2D *********************************** */

template <typename T, template <typename U> class Descriptor>
CountParticlesSelectiveFunctional2D<T, Descriptor>::CountParticlesSelectiveFunctional2D(
    util::SelectInt *tags_) :
    numParticlesId(this->getStatistics().subscribeIntSum()), tags(tags_)
{ }

template <typename T, template <typename U> class Descriptor>
CountParticlesSelectiveFunctional2D<T, Descriptor>::~CountParticlesSelectiveFunctional2D()
{
    delete tags;
}

template <typename T, template <typename U> class Descriptor>
CountParticlesSelectiveFunctional2D<T, Descriptor>::CountParticlesSelectiveFunctional2D(
    CountParticlesSelectiveFunctional2D<T, Descriptor> const &rhs) :
    numParticlesId(this->getStatistics().subscribeIntSum()), tags(rhs.tags->clone())
{ }

template <typename T, template <typename U> class Descriptor>
CountParticlesSelectiveFunctional2D<T, Descriptor>
    &CountParticlesSelectiveFunctional2D<T, Descriptor>::operator=(
        CountParticlesSelectiveFunctional2D<T, Descriptor> const &rhs)
{
    CountParticlesSelectiveFunctional2D<T, Descriptor>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void CountParticlesSelectiveFunctional2D<T, Descriptor>::swap(
    CountParticlesSelectiveFunctional2D<T, Descriptor> &rhs)
{
    std::swap(numParticlesId, rhs.numParticlesId);
    std::swap(tags, rhs.tags);
    PlainReductiveBoxProcessingFunctional2D::swap(rhs);
}

template <typename T, template <typename U> class Descriptor>
void CountParticlesSelectiveFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    std::vector<Particle2D<T, Descriptor> *> particles;
    particleField.findParticles(domain, particles);
    for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
        if ((*tags)(particles[iParticle]->getTag())) {
            this->getStatistics().gatherIntSum(numParticlesId, 1);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CountParticlesSelectiveFunctional2D<T, Descriptor>
    *CountParticlesSelectiveFunctional2D<T, Descriptor>::clone() const
{
    return new CountParticlesSelectiveFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CountParticlesSelectiveFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
plint CountParticlesSelectiveFunctional2D<T, Descriptor>::getNumParticles() const
{
    return this->getStatistics().getIntSum(numParticlesId);
}

/* ******** AverageParticleVelocityFunctional2D *********************************** */

template <typename T, template <typename U> class Descriptor>
AverageParticleVelocityFunctional2D<T, Descriptor>::AverageParticleVelocityFunctional2D() :
    averageVelocityId(Array<T, 2>(
        this->getStatistics().subscribeAverage(), this->getStatistics().subscribeAverage()))
{ }

template <typename T, template <typename U> class Descriptor>
void AverageParticleVelocityFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    std::vector<Particle2D<T, Descriptor> *> particles;
    particleField.findParticles(domain, particles);
    for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
        Array<T, 2> velocity;
        particles[iParticle]->getVelocity(velocity);
        this->getStatistics().gatherAverage(averageVelocityId[0], velocity[0]);
        this->getStatistics().gatherAverage(averageVelocityId[1], velocity[1]);
        this->getStatistics().incrementStats();
    }
}

template <typename T, template <typename U> class Descriptor>
AverageParticleVelocityFunctional2D<T, Descriptor>
    *AverageParticleVelocityFunctional2D<T, Descriptor>::clone() const
{
    return new AverageParticleVelocityFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void AverageParticleVelocityFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AverageParticleVelocityFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
Array<T, 2> AverageParticleVelocityFunctional2D<T, Descriptor>::getAverageParticleVelocity() const
{
    return Array<T, 2>(
        this->getStatistics().getAverage(averageVelocityId[0]),
        this->getStatistics().getAverage(averageVelocityId[1]));
}

/* ******** InjectParticlesFunctional2D *********************************** */

template <typename T, template <typename U> class Descriptor>
InjectParticlesFunctional2D<T, Descriptor>::InjectParticlesFunctional2D(
    std::vector<Particle2D<T, Descriptor> *> &particles_) :
    particles(particles_)
{ }

template <typename T, template <typename U> class Descriptor>
void InjectParticlesFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
        // The function addParticle tests automatically if the particle corresponds to the domain
        //   and needs to be added. In either case, the particle is consumed by the function.
        particleField.addParticle(domain, particles[iParticle]->clone());
    }
}

template <typename T, template <typename U> class Descriptor>
InjectParticlesFunctional2D<T, Descriptor> *InjectParticlesFunctional2D<T, Descriptor>::clone()
    const
{
    return new InjectParticlesFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InjectParticlesFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void InjectParticlesFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** InjectRandomPointParticlesFunctional2D *********************************** */

template <typename T, template <typename U> class Descriptor>
InjectRandomPointParticlesFunctional2D<T, Descriptor>::InjectRandomPointParticlesFunctional2D(
    plint tag_, T probabilityPerCell_) :
    tag(tag_), probabilityPerCell(probabilityPerCell_)
{ }

template <typename T, template <typename U> class Descriptor>
void InjectRandomPointParticlesFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T randNumber = (T)rand() / (T)RAND_MAX;
            if (randNumber < probabilityPerCell) {
                T randX = (T)rand() / (T)RAND_MAX - (T)1;
                T randY = (T)rand() / (T)RAND_MAX - (T)1;
                Array<T, 2> position(
                    particleField.getLocation().x + iX + randX,
                    particleField.getLocation().y + iY + randY);
                Array<T, 2> velocity;
                velocity.resetToZero();
                particleField.addParticle(
                    domain, new PointParticle2D<T, Descriptor>(tag, position, velocity));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
InjectRandomPointParticlesFunctional2D<T, Descriptor>
    *InjectRandomPointParticlesFunctional2D<T, Descriptor>::clone() const
{
    return new InjectRandomPointParticlesFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InjectRandomPointParticlesFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void InjectRandomPointParticlesFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** InjectRandomPointParticlesFunctionalPPRNG2D *********************************** */

template <typename T, template <typename U> class Descriptor>
InjectRandomPointParticlesFunctionalPPRNG2D<T, Descriptor>::
    InjectRandomPointParticlesFunctionalPPRNG2D(
        plint tag_, T probabilityPerCell_, Box2D const &boundingBox_, uint32_t const *seed_) :
    tag(tag_), probabilityPerCell(probabilityPerCell_), nY(boundingBox_.getNy()), seed(seed_)
{ }

template <typename T, template <typename U> class Descriptor>
void InjectRandomPointParticlesFunctionalPPRNG2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    Dot2D location = particleField.getLocation();
    sitmo::prng_engine eng(*seed);
    plint rng_index = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint globalX = nY * (iX + location.x);
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint globalY = iY + location.y + globalX;
            PLB_ASSERT(globalY >= rng_index);
            if (globalY > rng_index) {
                eng.discard(globalY - rng_index);
                rng_index = globalY;
            }
            T randNumber = (T)eng() / (T)sitmo::prng_engine::max();
            ++rng_index;
            if (randNumber < probabilityPerCell) {
                sitmo::prng_engine eng2(*seed + (uint32_t)globalY);
                T randX = (T)eng2() / (T)sitmo::prng_engine::max() - (T)0.5;
                T randY = (T)eng2() / (T)sitmo::prng_engine::max() - (T)0.5;
                Array<T, 2> position(location.x + iX + randX, location.y + iY + randY);
                Array<T, 2> velocity;
                velocity.resetToZero();
                particleField.addParticle(
                    domain, new PointParticle2D<T, Descriptor>(tag, position, velocity));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
InjectRandomPointParticlesFunctionalPPRNG2D<T, Descriptor>
    *InjectRandomPointParticlesFunctionalPPRNG2D<T, Descriptor>::clone() const
{
    return new InjectRandomPointParticlesFunctionalPPRNG2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InjectRandomPointParticlesFunctionalPPRNG2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void InjectRandomPointParticlesFunctionalPPRNG2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** AnalyticalInjectRandomParticlesFunctional2D *********************************** */

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>::
    AnalyticalInjectRandomParticlesFunctional2D(
        Particle2D<T, Descriptor> *particleTemplate_, T probabilityPerCell_,
        DomainFunctional functional_) :
    particleTemplate(particleTemplate_),
    probabilityPerCell(probabilityPerCell_),
    functional(functional_)
{ }

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>::
    AnalyticalInjectRandomParticlesFunctional2D(
        AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional> const &rhs) :
    particleTemplate(rhs.particleTemplate->clone()),
    probabilityPerCell(rhs.probabilityPerCell),
    functional(rhs.functional)
{ }

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>
    &AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>::operator=(
        AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional> const &rhs)
{
    AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
void AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>::swap(
    AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional> &rhs)
{
    std::swap(particleTemplate, rhs.particleTemplate);
    std::swap(probabilityPerCell, rhs.probabilityPerCell);
    std::swap(functional, rhs.functional);
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctional2D<
    T, Descriptor, DomainFunctional>::~AnalyticalInjectRandomParticlesFunctional2D()
{
    delete particleTemplate;
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
void AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>::
    processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T randNumber = (T)rand() / (T)RAND_MAX;
            if (randNumber < probabilityPerCell) {
                T randX = (T)rand() / (T)RAND_MAX - (T)1;
                T randY = (T)rand() / (T)RAND_MAX - (T)1;
                Array<T, 2> position(
                    particleField.getLocation().x + iX + randX,
                    particleField.getLocation().y + iY + randY);
                if (functional(position)) {
                    Particle2D<T, Descriptor> *newparticle = particleTemplate->clone();
                    newparticle->getPosition() = position;
                    particleField.addParticle(domain, newparticle);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>
    *AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>::clone() const
{
    return new AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>(*this);
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
void AnalyticalInjectRandomParticlesFunctional2D<T, Descriptor, DomainFunctional>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** AnalyticalInjectRandomParticlesFunctionalPPRNG2D *********************************** */

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>::
    AnalyticalInjectRandomParticlesFunctionalPPRNG2D(
        Particle2D<T, Descriptor> *particleTemplate_, T probabilityPerCell_,
        DomainFunctional functional_, Box2D const &boundingBox_, uint32_t const *seed_) :
    particleTemplate(particleTemplate_),
    probabilityPerCell(probabilityPerCell_),
    functional(functional_),
    nY(boundingBox_.getNy()),
    seed(seed_)
{ }

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>::
    AnalyticalInjectRandomParticlesFunctionalPPRNG2D(
        AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional> const
            &rhs) :
    particleTemplate(rhs.particleTemplate->clone()),
    probabilityPerCell(rhs.probabilityPerCell),
    functional(rhs.functional),
    nY(rhs.nY),
    seed(rhs.seed)
{ }

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>
    &AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>::operator=(
        AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional> const
            &rhs)
{
    AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>(rhs).swap(
        *this);
    return *this;
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
void AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>::swap(
    AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional> &rhs)
{
    std::swap(particleTemplate, rhs.particleTemplate);
    std::swap(probabilityPerCell, rhs.probabilityPerCell);
    std::swap(functional, rhs.functional);
    std::swap(nY, rhs.nY);
    std::swap(seed, rhs.seed);
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctionalPPRNG2D<
    T, Descriptor, DomainFunctional>::~AnalyticalInjectRandomParticlesFunctionalPPRNG2D()
{
    delete particleTemplate;
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
void AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>::
    processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    Dot2D location = particleField.getLocation();
    sitmo::prng_engine eng(*seed);
    plint rng_index = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint globalX = nY * (iX + location.x);
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint globalY = iY + location.y + globalX;
            PLB_ASSERT(globalY >= rng_index);
            if (globalY > rng_index) {
                eng.discard(globalY - rng_index);
                rng_index = globalY;
            }
            T randNumber = (T)eng() / (T)sitmo::prng_engine::max();
            ++rng_index;
            if (randNumber < probabilityPerCell) {
                sitmo::prng_engine eng2(*seed + (uint32_t)globalY);
                T randX = (T)eng2() / (T)sitmo::prng_engine::max() - (T)0.5;
                T randY = (T)eng2() / (T)sitmo::prng_engine::max() - (T)0.5;
                Array<T, 2> position(location.x + iX + randX, location.y + iY + randY);
                if (functional(position)) {
                    Particle2D<T, Descriptor> *newparticle = particleTemplate->clone();
                    newparticle->getPosition() = position;
                    particleField.addParticle(domain, newparticle);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional> *
    AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>::clone() const
{
    return new AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>(
        *this);
}

template <typename T, template <typename U> class Descriptor, class DomainFunctional>
void AnalyticalInjectRandomParticlesFunctionalPPRNG2D<T, Descriptor, DomainFunctional>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** AbsorbParticlesFunctional2D *********************************** */

template <typename T, template <typename U> class Descriptor>
void AbsorbParticlesFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    particleField.removeParticles(domain);
}

template <typename T, template <typename U> class Descriptor>
AbsorbParticlesFunctional2D<T, Descriptor> *AbsorbParticlesFunctional2D<T, Descriptor>::clone()
    const
{
    return new AbsorbParticlesFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AbsorbParticlesFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void AbsorbParticlesFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** FluidToParticleCoupling2D *********************************** */

template <typename T, template <typename U> class Descriptor>
FluidToParticleCoupling2D<T, Descriptor>::FluidToParticleCoupling2D(T scaling_) : scaling(scaling_)
{ }

template <typename T, template <typename U> class Descriptor>
void FluidToParticleCoupling2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    BlockLattice2D<T, Descriptor> &fluid =
        *dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[1]);
    particleField.fluidToParticleCoupling(domain, fluid, scaling);
}

template <typename T, template <typename U> class Descriptor>
FluidToParticleCoupling2D<T, Descriptor> *FluidToParticleCoupling2D<T, Descriptor>::clone() const
{
    return new FluidToParticleCoupling2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT FluidToParticleCoupling2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;  // Important: must be bulk-only,
    // because non-local access to fluid is made.
}

template <typename T, template <typename U> class Descriptor>
void FluidToParticleCoupling2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;           // Fluid.
}

/* ******** VelocityToParticleCoupling2D *********************************** */

template <typename T, template <typename U> class Descriptor>
VelocityToParticleCoupling2D<T, Descriptor>::VelocityToParticleCoupling2D(T scaling_) :
    scaling(scaling_)
{ }

template <typename T, template <typename U> class Descriptor>
void VelocityToParticleCoupling2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField2D<T, Descriptor> *particleFieldPtr =
        dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(particleFieldPtr);
    ParticleField2D<T, Descriptor> &particleField = *particleFieldPtr;

    TensorField2D<T, 2> *velocityPtr = dynamic_cast<TensorField2D<T, 2> *>(blocks[1]);
    PLB_ASSERT(velocityPtr);
    TensorField2D<T, 2> &velocity = *velocityPtr;

    particleField.velocityToParticleCoupling(domain, velocity, scaling);
}

template <typename T, template <typename U> class Descriptor>
VelocityToParticleCoupling2D<T, Descriptor> *VelocityToParticleCoupling2D<T, Descriptor>::clone()
    const
{
    return new VelocityToParticleCoupling2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void VelocityToParticleCoupling2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;           // TensorField velocity.
}

/* ******** N_VelocityToParticleCoupling2D *********************************** */

template <typename T, template <typename U> class Descriptor>
N_VelocityToParticleCoupling2D<T, Descriptor>::N_VelocityToParticleCoupling2D(T scaling_) :
    scaling(scaling_)
{ }

template <typename T, template <typename U> class Descriptor>
void N_VelocityToParticleCoupling2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField2D<T, Descriptor> *particleFieldPtr =
        dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(particleFieldPtr);
    ParticleField2D<T, Descriptor> &particleField = *particleFieldPtr;

    NTensorField2D<T> *velocityPtr = dynamic_cast<NTensorField2D<T> *>(blocks[1]);
    PLB_ASSERT(velocityPtr);
    NTensorField2D<T> &velocity = *velocityPtr;
    PLB_ASSERT(velocity.getNdim() == 2);

    particleField.velocityToParticleCoupling(domain, velocity, scaling);
}

template <typename T, template <typename U> class Descriptor>
N_VelocityToParticleCoupling2D<T, Descriptor>
    *N_VelocityToParticleCoupling2D<T, Descriptor>::clone() const
{
    return new N_VelocityToParticleCoupling2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void N_VelocityToParticleCoupling2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;           // TensorField velocity.
}

/* ******** RhoBarJtoParticleCoupling2D *********************************** */

template <typename T, template <typename U> class Descriptor>
RhoBarJtoParticleCoupling2D<T, Descriptor>::RhoBarJtoParticleCoupling2D(bool velIsJ_, T scaling_) :
    velIsJ(velIsJ_), scaling(scaling_)
{ }

template <typename T, template <typename U> class Descriptor>
void RhoBarJtoParticleCoupling2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField2D<T, Descriptor> *particleFieldPtr =
        dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(particleFieldPtr);
    ParticleField2D<T, Descriptor> &particleField = *particleFieldPtr;

    NTensorField2D<T> *rhoBarJptr = dynamic_cast<NTensorField2D<T> *>(blocks[1]);
    PLB_ASSERT(rhoBarJptr);
    NTensorField2D<T> &rhoBarJ = *rhoBarJptr;

    particleField.rhoBarJtoParticleCoupling(domain, rhoBarJ, velIsJ, scaling);
}

template <typename T, template <typename U> class Descriptor>
RhoBarJtoParticleCoupling2D<T, Descriptor> *RhoBarJtoParticleCoupling2D<T, Descriptor>::clone()
    const
{
    return new RhoBarJtoParticleCoupling2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void RhoBarJtoParticleCoupling2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;           // RhoBarJ;
}

/* ******** AdvanceParticlesFunctional2D *********************************** */

template <typename T, template <typename U> class Descriptor>
AdvanceParticlesFunctional2D<T, Descriptor>::AdvanceParticlesFunctional2D(T cutOffValue_) :
    cutOffValue(cutOffValue_)
{ }

template <typename T, template <typename U> class Descriptor>
void AdvanceParticlesFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    particleField.advanceParticles(domain, cutOffValue);
}

template <typename T, template <typename U> class Descriptor>
AdvanceParticlesFunctional2D<T, Descriptor> *AdvanceParticlesFunctional2D<T, Descriptor>::clone()
    const
{
    return new AdvanceParticlesFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT AdvanceParticlesFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;  // Important: access envelope as well,
                                          // because particles are streamed from the
                                          // envelope into the bulk.
}

template <typename T, template <typename U> class Descriptor>
void AdvanceParticlesFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** AdvanceParticlesFunctional2D *********************************** */

template <typename T, template <typename U> class Descriptor>
AdvanceParticlesEveryWhereFunctional2D<T, Descriptor>::AdvanceParticlesEveryWhereFunctional2D(
    T cutOffValue_) :
    cutOffValue(cutOffValue_)
{ }

template <typename T, template <typename U> class Descriptor>
void AdvanceParticlesEveryWhereFunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    // Important: acts on the full domain, not only the provided one.
    particleField.advanceParticles(particleField.getBoundingBox(), cutOffValue);
}

template <typename T, template <typename U> class Descriptor>
AdvanceParticlesEveryWhereFunctional2D<T, Descriptor>
    *AdvanceParticlesEveryWhereFunctional2D<T, Descriptor>::clone() const
{
    return new AdvanceParticlesEveryWhereFunctional2D<T, Descriptor>(*this);
}

// template<typename T, template<typename U> class Descriptor>
// BlockDomain::DomainT AdvanceParticlesEveryWhereFunctional2D<T,Descriptor>::appliesTo() const {
//     return BlockDomain::bulkAndEnvelope;  // Important: access envelope as well,
//                                           // because particles are streamed from the
//                                           // envelope into the bulk.}
//
// }

template <typename T, template <typename U> class Descriptor>
void AdvanceParticlesEveryWhereFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** VerletUpdateVelocity2D *********************************** */

template <typename T, template <typename U> class Descriptor>
void VerletUpdateVelocity2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);

    std::vector<Particle2D<T, Descriptor> *> found;
    particleField.findParticles(domain, found);

    for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
        Particle2D<T, Descriptor> *nonTypeParticle = found[iParticle];
        VerletParticle2D<T, Descriptor> *particle =
            dynamic_cast<VerletParticle2D<T, Descriptor> *>(nonTypeParticle);
        PLB_ASSERT(particle);

        Array<T, 2> a(particle->get_a());
        particle->set_v(particle->get_vHalfTime() + (T)0.5 * a);
        PLB_ASSERT(norm(particle->get_v()) < 1.);
    }
}

template <typename T, template <typename U> class Descriptor>
VerletUpdateVelocity2D<T, Descriptor> *VerletUpdateVelocity2D<T, Descriptor>::clone() const
{
    return new VerletUpdateVelocity2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void VerletUpdateVelocity2D<T, Descriptor>::getModificationPattern(
    std::vector<bool> &isWritten) const
{
    isWritten[0] = true;  // Particle field.
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT VerletUpdateVelocity2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void VerletUpdateVelocity2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}

/* ******** CountAndAccumulateParticles2D *********************************** */

template <typename T, template <typename U> class Descriptor>
void CountAndAccumulateParticles2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    ScalarField2D<plint> &numParticlefield = *dynamic_cast<ScalarField2D<plint> *>(blocks[1]);
    Dot2D offset = computeRelativeDisplacement(particleField, numParticlefield);
    std::vector<Particle2D<T, Descriptor> *> particles;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            particleField.findParticles(Box2D(iX, iX, iY, iY), particles);
            numParticlefield.get(iX + offset.x, iY + offset.y) += particles.size();
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CountAndAccumulateParticles2D<T, Descriptor> *CountAndAccumulateParticles2D<T, Descriptor>::clone()
    const
{
    return new CountAndAccumulateParticles2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT CountAndAccumulateParticles2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void CountAndAccumulateParticles2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Particle field.
    modified[1] = modif::staticVariables;  // Scalar field.
}

/* ******** CountAndAccumulateTaggedParticles2D *********************************** */

template <typename T, template <typename U> class Descriptor>
CountAndAccumulateTaggedParticles2D<T, Descriptor>::CountAndAccumulateTaggedParticles2D(
    plint tag_) :
    tag(tag_)
{ }

template <typename T, template <typename U> class Descriptor>
void CountAndAccumulateTaggedParticles2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    ScalarField2D<plint> &numParticlefield = *dynamic_cast<ScalarField2D<plint> *>(blocks[1]);
    Dot2D offset = computeRelativeDisplacement(particleField, numParticlefield);
    std::vector<Particle2D<T, Descriptor> *> particles;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            particleField.findParticles(Box2D(iX, iX, iY, iY), particles);
            for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
                if (particles[iParticle]->getTag() == tag) {
                    ++numParticlefield.get(iX + offset.x, iY + offset.y);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CountAndAccumulateTaggedParticles2D<T, Descriptor>
    *CountAndAccumulateTaggedParticles2D<T, Descriptor>::clone() const
{
    return new CountAndAccumulateTaggedParticles2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT CountAndAccumulateTaggedParticles2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void CountAndAccumulateTaggedParticles2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Particle field.
    modified[1] = modif::staticVariables;  // Scalar field.
}

/* ******** CountTaggedParticles2D *********************************** */

template <typename T, template <typename U> class Descriptor>
CountTaggedParticles2D<T, Descriptor>::CountTaggedParticles2D(util::SelectInt *tags_) : tags(tags_)
{ }

template <typename T, template <typename U> class Descriptor>
CountTaggedParticles2D<T, Descriptor>::~CountTaggedParticles2D()
{
    delete tags;
}

template <typename T, template <typename U> class Descriptor>
CountTaggedParticles2D<T, Descriptor>::CountTaggedParticles2D(
    CountTaggedParticles2D<T, Descriptor> const &rhs) :
    tags(rhs.tags->clone())
{ }

template <typename T, template <typename U> class Descriptor>
CountTaggedParticles2D<T, Descriptor> &CountTaggedParticles2D<T, Descriptor>::operator=(
    CountTaggedParticles2D<T, Descriptor> const &rhs)
{
    CountTaggedParticles2D<T, Descriptor>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void CountTaggedParticles2D<T, Descriptor>::swap(CountTaggedParticles2D<T, Descriptor> &rhs)
{
    std::swap(tags, rhs.tags);
}

template <typename T, template <typename U> class Descriptor>
void CountTaggedParticles2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField2D<T, Descriptor> &particleField =
        *dynamic_cast<ParticleField2D<T, Descriptor> *>(blocks[0]);
    ScalarField2D<plint> &numParticlefield = *dynamic_cast<ScalarField2D<plint> *>(blocks[1]);
    Dot2D offset = computeRelativeDisplacement(particleField, numParticlefield);
    std::vector<Particle2D<T, Descriptor> *> particles;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            particleField.findParticles(Box2D(iX, iX, iY, iY), particles);
            numParticlefield.get(iX + offset.x, iY + offset.y) = 0;
            for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
                if ((*tags)(particles[iParticle]->getTag())) {
                    ++numParticlefield.get(iX + offset.x, iY + offset.y);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CountTaggedParticles2D<T, Descriptor> *CountTaggedParticles2D<T, Descriptor>::clone() const
{
    return new CountTaggedParticles2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CountTaggedParticles2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Particle field.
    modified[1] = modif::staticVariables;  // Scalar field.
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T_, template <typename U_> class Descriptor_> class ParticleFieldT>
plint countParticles(
    MultiParticleField2D<ParticleFieldT<T, Descriptor> > &particles, Box2D const &domain)
{
    std::vector<MultiBlock2D *> particleArg;
    particleArg.push_back(&particles);

    CountParticlesFunctional2D<T, Descriptor> functional;
    applyProcessingFunctional(functional, domain, particleArg);
    return functional.getNumParticles();
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T_, template <typename U_> class Descriptor_> class ParticleFieldT>
plint countParticles(
    MultiParticleField2D<ParticleFieldT<T, Descriptor> > &particles, Box2D const &domain,
    util::SelectInt *tags)
{
    std::vector<MultiBlock2D *> particleArg;
    particleArg.push_back(&particles);

    CountParticlesSelectiveFunctional2D<T, Descriptor> functional(tags);
    applyProcessingFunctional(functional, domain, particleArg);
    return functional.getNumParticles();
}

template <typename T, template <typename U> class Descriptor>
void injectParticles(
    std::vector<Particle2D<T, Descriptor> *> &injectedParticles,
    MultiParticleField2D<DenseParticleField2D<T, Descriptor> > &particles, Box2D domain)
{
    std::vector<MultiBlock2D *> particleArg;
    particleArg.push_back(&particles);
    applyProcessingFunctional(
        new InjectParticlesFunctional2D<T, Descriptor>(injectedParticles), domain, particleArg);
}

}  // namespace plb

#endif  // PARTICLE_PROCESSING_FUNCTIONAL_2D_HH
