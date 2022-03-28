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

#ifndef PARTICLE_FIELD_2D_HH
#define PARTICLE_FIELD_2D_HH

#include <utility>

#include "core/globalDefs.h"
#include "particles/particleField2D.h"

namespace plb {

/* *************** class ParticleField2D ************************************ */

template <typename T, template <typename U> class Descriptor>
ParticleField2D<T, Descriptor>::ParticleField2D(plint nx, plint ny) : AtomicBlock2D(nx, ny)
{ }

template <typename T, template <typename U> class Descriptor>
bool ParticleField2D<T, Descriptor>::isContained(Array<T, 2> const &particlePos, Box2D box) const
{
    Dot2D const &location = this->getLocation();
    T x = particlePos[0] - location.x;
    T y = particlePos[1] - location.y;

    return (x > (T)box.x0 - (T)0.5) && (x <= (T)box.x1 + (T)0.5) && (y > (T)box.y0 - (T)0.5)
           && (y <= (T)box.y1 + (T)0.5);
}

template <typename T, template <typename U> class Descriptor>
plint ParticleField2D<T, Descriptor>::nearestCell(T pos)
{
    T afterComma = pos - floor(pos);
    if (pos >= (T)0.) {
        if (afterComma > (T)0.5) {
            return (plint)(pos + (T)0.75);
        } else {
            return (plint)(pos);
        }
    } else {
        if (afterComma > (T)0.5) {
            return (plint)(pos);
        } else {
            return (plint)(pos - (T)0.75);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void ParticleField2D<T, Descriptor>::computeGridPosition(
    Array<T, 2> const &position, plint &iX, plint &iY) const
{
    Dot2D const &location = this->getLocation();
    iX = nearestCell(position[0]) - location.x;
    iY = nearestCell(position[1]) - location.y;
}

/* *************** class DenseParticleDataTransfer2D ************************ */

template <typename T, template <typename U> class Descriptor>
DenseParticleDataTransfer2D<T, Descriptor>::DenseParticleDataTransfer2D(
    DenseParticleField2D<T, Descriptor> &particleField_) :
    particleField(particleField_)
{ }

template <typename T, template <typename U> class Descriptor>
plint DenseParticleDataTransfer2D<T, Descriptor>::staticCellSize() const
{
    return 0;  // Particle containers have only dynamic data.
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer2D<T, Descriptor>::send(
    Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    buffer.clear();
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the send procedure
    //   is run whenever kind is one of the dynamic types.
    if ((kind == modif::dynamicVariables) || (kind == modif::allVariables)
        || (kind == modif::dataStructure))
    {
        std::vector<Particle2D<T, Descriptor> *> foundParticles;
        particleField.findParticles(domain, foundParticles);
        for (pluint iParticle = 0; iParticle < foundParticles.size(); ++iParticle) {
            // The serialize function automatically reallocates memory for buffer.
            serialize(*foundParticles[iParticle], buffer);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer2D<T, Descriptor>::receive(
    Box2D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(contained(domain, particleField.getBoundingBox()));
    // Clear the existing data before introducing the new data.
    particleField.removeParticles(domain);
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the receive procedure
    //   is run whenever kind is one of the dynamic types.
    if ((kind == modif::dynamicVariables) || (kind == modif::allVariables)
        || (kind == modif::dataStructure))
    {
        pluint posInBuffer = 0;
        while (posInBuffer < buffer.size()) {
            // 1. Generate dynamics object, and unserialize dynamic data.
            HierarchicUnserializer unserializer(buffer, posInBuffer);
            Particle2D<T, Descriptor> *newParticle =
                meta::particleRegistration2D<T, Descriptor>().generate(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            particleField.addParticle(domain, newParticle);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer2D<T, Descriptor>::receive(
    Box2D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot2D absoluteOffset)
{
    if (absoluteOffset.x == 0 && absoluteOffset.y == 0) {
        receive(domain, buffer, kind);
        return;
    }
    PLB_PRECONDITION(contained(domain, particleField.getBoundingBox()));
    Array<T, 2> realAbsoluteOffset((T)absoluteOffset.x, (T)absoluteOffset.y);
    // Clear the existing data before introducing the new data.
    particleField.removeParticles(domain);
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the receive procedure
    //   is run whenever kind is one of the dynamic types.
    if ((kind == modif::dynamicVariables) || (kind == modif::allVariables)
        || (kind == modif::dataStructure))
    {
        pluint posInBuffer = 0;
        while (posInBuffer < buffer.size()) {
            // 1. Generate dynamics object, and unserialize dynamic data.
            HierarchicUnserializer unserializer(buffer, posInBuffer);
            Particle2D<T, Descriptor> *newParticle =
                meta::particleRegistration2D<T, Descriptor>().generate(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            newParticle->getPosition() += realAbsoluteOffset;

            particleField.addParticle(domain, newParticle);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer2D<T, Descriptor>::attribute(
    Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind)
{
    Box2D fromDomain(toDomain.shift(deltaX, deltaY));
    std::vector<char> buffer;
    DenseParticleField2D<T, Descriptor> const &fromParticleField =
        dynamic_cast<DenseParticleField2D<T, Descriptor> const &>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind);
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer2D<T, Descriptor>::attribute(
    Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind,
    Dot2D absoluteOffset)
{
    Box2D fromDomain(toDomain.shift(deltaX, deltaY));
    std::vector<char> buffer;
    DenseParticleField2D<T, Descriptor> const &fromParticleField =
        dynamic_cast<DenseParticleField2D<T, Descriptor> const &>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind, absoluteOffset);
}

/* *************** class DenseParticleField2D ********************** */

template <typename T, template <typename U> class Descriptor>
DenseParticleField2D<T, Descriptor>::DenseParticleField2D(plint nx, plint ny) :
    ParticleField2D<T, Descriptor>(nx, ny), particleGrid(nx, ny), dataTransfer(*this)
{ }

template <typename T, template <typename U> class Descriptor>
DenseParticleField2D<T, Descriptor>::~DenseParticleField2D()
{
    for (plint iX = 0; iX < particleGrid.getNx(); ++iX) {
        for (plint iY = 0; iY < particleGrid.getNy(); ++iY) {
            for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY).size(); ++iParticle) {
                delete particleGrid.get(iX, iY)[iParticle];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
DenseParticleField2D<T, Descriptor>::DenseParticleField2D(DenseParticleField2D const &rhs) :
    ParticleField2D<T, Descriptor>(rhs),
    particleGrid(rhs.particleGrid.getNx(), rhs.particleGrid.getNy()),
    dataTransfer(*this)

{
    for (plint iX = 0; iX < particleGrid.getNx(); ++iX) {
        for (plint iY = 0; iY < particleGrid.getNy(); ++iY) {
            particleGrid.get(iX, iY).resize(rhs.particleGrid.get(iX, iY).size());
            for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY).size(); ++iParticle) {
                particleGrid.get(iX, iY)[iParticle] =
                    rhs.particleGrid.get(iX, iY)[iParticle]->clone();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
DenseParticleField2D<T, Descriptor> &DenseParticleField2D<T, Descriptor>::operator=(
    DenseParticleField2D<T, Descriptor> const &rhs)
{
    DenseParticleField2D<T, Descriptor>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
DenseParticleField2D<T, Descriptor> *DenseParticleField2D<T, Descriptor>::clone() const
{
    return new DenseParticleField2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::swap(DenseParticleField2D<T, Descriptor> &rhs)
{
    ParticleField2D<T, Descriptor>::swap(rhs);
    particleGrid.swap(rhs.particleGrid);
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::addParticle(
    Box2D domain, Particle2D<T, Descriptor> *particle)
{
    plint iX, iY;
    this->computeGridPosition(particle->getPosition(), iX, iY);
    Box2D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)
        && contained(iX, iY, finalDomain)) {
        particleGrid.get(iX, iY).push_back(particle);
    } else {
        delete particle;
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::removeParticles(Box2D domain)
{
    Box2D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY).size(); ++iParticle)
                {
                    delete particleGrid.get(iX, iY)[iParticle];
                }
                particleGrid.get(iX, iY).clear();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::removeParticles(Box2D domain, plint tag)
{
    Box2D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                typename std::vector<Particle2D<T, Descriptor> *>::iterator it =
                    particleGrid.get(iX, iY).begin();
                for (; it != particleGrid.get(iX, iY).end();) {
                    if ((*it)->getTag() == tag) {
                        delete *it;
                        it = particleGrid.get(iX, iY).erase(it);
                    } else {
                        ++it;
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::findParticles(
    Box2D domain, std::vector<Particle2D<T, Descriptor> *> &found)
{
    found.clear();
    PLB_ASSERT(contained(domain, particleGrid.getBoundingBox()));
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY).size(); ++iParticle) {
                found.push_back(particleGrid.get(iX, iY)[iParticle]);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::findParticles(
    Box2D domain, std::vector<Particle2D<T, Descriptor> const *> &found) const
{
    found.clear();
    PLB_ASSERT(contained(domain, particleGrid.getBoundingBox()));
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY).size(); ++iParticle) {
                found.push_back(particleGrid.get(iX, iY)[iParticle]);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::velocityToParticleCoupling(
    Box2D domain, TensorField2D<T, 2> &velocityField, T scaling)
{
    Box2D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY).size(); ++iParticle)
                {
                    particleGrid.get(iX, iY)[iParticle]->velocityToParticle(velocityField, scaling);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::velocityToParticleCoupling(
    Box2D domain, NTensorField2D<T> &velocityField, T scaling)
{
    Box2D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY).size(); ++iParticle)
                {
                    particleGrid.get(iX, iY)[iParticle]->velocityToParticle(velocityField, scaling);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::rhoBarJtoParticleCoupling(
    Box2D domain, NTensorField2D<T> &rhoBarJfield, bool velIsJ, T scaling)
{
    Box2D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY).size(); ++iParticle)
                {
                    particleGrid.get(iX, iY)[iParticle]->rhoBarJtoParticle(
                        rhoBarJfield, velIsJ, scaling);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::fluidToParticleCoupling(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, T scaling)
{
    Box2D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY).size(); ++iParticle)
                {
                    particleGrid.get(iX, iY)[iParticle]->fluidToParticle(lattice, scaling);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField2D<T, Descriptor>::advanceParticles(Box2D domain, T cutOffValue)
{
    Box2D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        std::vector<std::pair<Dot2D, Particle2D<T, Descriptor> *> > nextCellParticles;
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                std::vector<Particle2D<T, Descriptor> *> &particles = particleGrid.get(iX, iY);
                std::vector<Particle2D<T, Descriptor> *> newLocalParticles;
                for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
                    Particle2D<T, Descriptor> *particle = particles[iParticle];
                    Array<T, 2> oldPos(particle->getPosition());
                    particle->advance();
                    if (cutOffValue >= T()
                        && normSqr(oldPos - particle->getPosition()) < cutOffValue) {
                        delete particle;
                    } else {
                        plint newX, newY;
                        this->computeGridPosition(particle->getPosition(), newX, newY);
                        if (newX == iX && newY == iY) {
                            newLocalParticles.push_back(particle);
                        } else {
                            if (contained(newX, newY, finalDomain)) {
                                nextCellParticles.push_back(
                                    std::make_pair(Dot2D(newX, newY), particle));
                            } else {
                                delete particle;
                            }
                        }
                    }
                }
                newLocalParticles.swap(particles);
            }
        }
        for (pluint i = 0; i < nextCellParticles.size(); ++i) {
            plint newX = nextCellParticles[i].first.x;
            plint newY = nextCellParticles[i].first.y;
            particleGrid.get(newX, newY).push_back(nextCellParticles[i].second);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
DenseParticleDataTransfer2D<T, Descriptor> &DenseParticleField2D<T, Descriptor>::getDataTransfer()
{
    return dataTransfer;
}

template <typename T, template <typename U> class Descriptor>
DenseParticleDataTransfer2D<T, Descriptor> const &
    DenseParticleField2D<T, Descriptor>::getDataTransfer() const
{
    return dataTransfer;
}

template <typename T, template <typename U> class Descriptor>
std::string DenseParticleField2D<T, Descriptor>::getBlockName()
{
    return std::string("DenseParticleField2D");
}

template <typename T, template <typename U> class Descriptor>
std::string DenseParticleField2D<T, Descriptor>::basicType()
{
    return std::string(NativeType<T>::getName());
}

template <typename T, template <typename U> class Descriptor>
std::string DenseParticleField2D<T, Descriptor>::descriptorType()
{
    return std::string(Descriptor<T>::name);
}

/* *************** class LightParticleDataTransfer2D ************************ */

template <typename T, template <typename U> class Descriptor>
LightParticleDataTransfer2D<T, Descriptor>::LightParticleDataTransfer2D(
    LightParticleField2D<T, Descriptor> &particleField_) :
    particleField(particleField_)
{ }

template <typename T, template <typename U> class Descriptor>
plint LightParticleDataTransfer2D<T, Descriptor>::staticCellSize() const
{
    return 0;  // Particle containers have only dynamic data.
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer2D<T, Descriptor>::send(
    Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    buffer.clear();
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the send procedure
    //   is run whenever kind is one of the dynamic types.
    if ((kind == modif::dynamicVariables) || (kind == modif::allVariables)
        || (kind == modif::dataStructure))
    {
        std::vector<Particle2D<T, Descriptor> *> foundParticles;
        particleField.findParticles(domain, foundParticles);
        for (pluint iParticle = 0; iParticle < foundParticles.size(); ++iParticle) {
            // The serialize function automatically reallocates memory for buffer.
            serialize(*foundParticles[iParticle], buffer);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer2D<T, Descriptor>::receive(
    Box2D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(contained(domain, particleField.getBoundingBox()));
    // Clear the existing data before introducing the new data.
    particleField.removeParticles(domain);
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the receive procedure
    //   is run whenever kind is one of the dynamic types.
    if ((kind == modif::dynamicVariables) || (kind == modif::allVariables)
        || (kind == modif::dataStructure))
    {
        pluint posInBuffer = 0;
        while (posInBuffer < buffer.size()) {
            // 1. Generate dynamics object, and unserialize dynamic data.
            HierarchicUnserializer unserializer(buffer, posInBuffer);
            Particle2D<T, Descriptor> *newParticle =
                meta::particleRegistration2D<T, Descriptor>().generate(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            particleField.addParticle(domain, newParticle);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer2D<T, Descriptor>::receive(
    Box2D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot2D absoluteOffset)
{
    if (absoluteOffset.x == 0 && absoluteOffset.y == 0) {
        receive(domain, buffer, kind);
        return;
    }
    PLB_PRECONDITION(contained(domain, particleField.getBoundingBox()));
    Array<T, 3> realAbsoluteOffset((T)absoluteOffset.x, (T)absoluteOffset.y);
    // Clear the existing data before introducing the new data.
    particleField.removeParticles(domain);
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the receive procedure
    //   is run whenever kind is one of the dynamic types.
    if ((kind == modif::dynamicVariables) || (kind == modif::allVariables)
        || (kind == modif::dataStructure))
    {
        pluint posInBuffer = 0;
        while (posInBuffer < buffer.size()) {
            // 1. Generate dynamics object, and unserialize dynamic data.
            HierarchicUnserializer unserializer(buffer, posInBuffer);
            Particle2D<T, Descriptor> *newParticle =
                meta::particleRegistration2D<T, Descriptor>().generate(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            newParticle->getPosition() += realAbsoluteOffset;

            particleField.addParticle(domain, newParticle);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer2D<T, Descriptor>::attribute(
    Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind)
{
    Box2D fromDomain(toDomain.shift(deltaX, deltaY));
    std::vector<char> buffer;
    LightParticleField2D<T, Descriptor> const &fromParticleField =
        dynamic_cast<LightParticleField2D<T, Descriptor> const &>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind);
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer2D<T, Descriptor>::attribute(
    Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind,
    Dot2D absoluteOffset)
{
    Box2D fromDomain(toDomain.shift(deltaX, deltaY));
    std::vector<char> buffer;
    LightParticleField2D<T, Descriptor> const &fromParticleField =
        dynamic_cast<LightParticleField2D<T, Descriptor> const &>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind, absoluteOffset);
}

/* *************** class LightParticleField2D ********************** */

template <typename T, template <typename U> class Descriptor>
LightParticleField2D<T, Descriptor>::LightParticleField2D(plint nx, plint ny) :
    ParticleField2D<T, Descriptor>(nx, ny), dataTransfer(*this)
{ }

template <typename T, template <typename U> class Descriptor>
LightParticleField2D<T, Descriptor>::~LightParticleField2D()
{
    for (pluint i = 0; i < particles.size(); ++i) {
        delete particles[i];
    }
}

template <typename T, template <typename U> class Descriptor>
LightParticleField2D<T, Descriptor>::LightParticleField2D(LightParticleField2D const &rhs) :
    ParticleField2D<T, Descriptor>(rhs), dataTransfer(*this)
{
    for (pluint i = 0; i < rhs.particles.size(); ++i) {
        particles.push_back(rhs.particles[i].clone());
    }
}

template <typename T, template <typename U> class Descriptor>
LightParticleField2D<T, Descriptor> &LightParticleField2D<T, Descriptor>::operator=(
    LightParticleField2D<T, Descriptor> const &rhs)
{
    LightParticleField2D<T, Descriptor>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
LightParticleField2D<T, Descriptor> *LightParticleField2D<T, Descriptor>::clone() const
{
    return new LightParticleField2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::swap(LightParticleField2D<T, Descriptor> &rhs)
{
    ParticleField2D<T, Descriptor>::swap(rhs);
    particles.swap(rhs.particles);
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::addParticle(
    Box2D domain, Particle2D<T, Descriptor> *particle)
{
    Box2D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)
        && this->isContained(particle->getPosition(), finalDomain))
    {
        particles.push_back(particle);
    } else {
        delete particle;
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::removeParticles(Box2D domain)
{
    std::vector<Particle2D<T, Descriptor> *> remainingParticles;
    Box2D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            if (this->isContained(particles[i]->getPosition(), finalDomain)) {
                delete particles[i];
            } else {
                remainingParticles.push_back(particles[i]);
            }
        }
    }
    remainingParticles.swap(particles);
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::removeParticles(Box2D domain, plint tag)
{
    std::vector<Particle2D<T, Descriptor> *> remainingParticles;
    Box2D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            if (this->isContained(particles[i]->getPosition(), finalDomain)
                && particles[i]->getTag() == tag) {
                delete particles[i];
            } else {
                remainingParticles.push_back(particles[i]);
            }
        }
    }
    remainingParticles.swap(particles);
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::findParticles(
    Box2D domain, std::vector<Particle2D<T, Descriptor> *> &found)
{
    found.clear();
    PLB_ASSERT(contained(domain, this->getBoundingBox()));
    for (pluint i = 0; i < particles.size(); ++i) {
        if (this->isContained(particles[i]->getPosition(), domain)) {
            found.push_back(particles[i]);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::findParticles(
    Box2D domain, std::vector<Particle2D<T, Descriptor> const *> &found) const
{
    found.clear();
    PLB_ASSERT(contained(domain, this->getBoundingBox()));
    for (pluint i = 0; i < particles.size(); ++i) {
        if (this->isContained(particles[i]->getPosition(), domain)) {
            found.push_back(particles[i]);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::velocityToParticleCoupling(
    Box2D domain, TensorField2D<T, 2> &velocityField, T scaling)
{
    Box2D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            if (this->isContained(particles[i]->getPosition(), finalDomain)) {
                particles[i]->velocityToParticle(velocityField, scaling);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::velocityToParticleCoupling(
    Box2D domain, NTensorField2D<T> &velocityField, T scaling)
{
    Box2D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            if (this->isContained(particles[i]->getPosition(), finalDomain)) {
                particles[i]->velocityToParticle(velocityField, scaling);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::rhoBarJtoParticleCoupling(
    Box2D domain, NTensorField2D<T> &rhoBarJfield, bool velIsJ, T scaling)
{
    Box2D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            if (this->isContained(particles[i]->getPosition(), finalDomain)) {
                particles[i]->rhoBarJtoParticle(rhoBarJfield, velIsJ, scaling);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::fluidToParticleCoupling(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, T scaling)
{
    Box2D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            if (this->isContained(particles[i]->getPosition(), finalDomain)) {
                particles[i]->fluidToParticle(lattice, scaling);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField2D<T, Descriptor>::advanceParticles(Box2D domain, T cutOffValue)
{
    std::vector<Particle2D<T, Descriptor> *> remainingParticles;
    Box2D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            Particle2D<T, Descriptor> *particle = particles[i];
            if (this->isContained(particle->getPosition(), finalDomain)) {
                Array<T, 3> oldPos(particle->getPosition());
                particle->advance();
                if ((cutOffValue >= T() && normSqr(oldPos - particle->getPosition()) < cutOffValue)
                    || (!this->isContained(particle->getPosition(), this->getBoundingBox())))
                {
                    delete particle;
                } else {
                    remainingParticles.push_back(particle);
                }
            }
        }
    }
    particles.swap(remainingParticles);
}

template <typename T, template <typename U> class Descriptor>
LightParticleDataTransfer2D<T, Descriptor> &LightParticleField2D<T, Descriptor>::getDataTransfer()
{
    return dataTransfer;
}

template <typename T, template <typename U> class Descriptor>
LightParticleDataTransfer2D<T, Descriptor> const &
    LightParticleField2D<T, Descriptor>::getDataTransfer() const
{
    return dataTransfer;
}

template <typename T, template <typename U> class Descriptor>
std::string LightParticleField2D<T, Descriptor>::getBlockName()
{
    return std::string("LightParticleField2D");
}

template <typename T, template <typename U> class Descriptor>
std::string LightParticleField2D<T, Descriptor>::basicType()
{
    return std::string(NativeType<T>::getName());
}

template <typename T, template <typename U> class Descriptor>
std::string LightParticleField2D<T, Descriptor>::descriptorType()
{
    return std::string(Descriptor<T>::name);
}

}  // namespace plb

#endif  // PARTICLE_FIELD_2D_HH
