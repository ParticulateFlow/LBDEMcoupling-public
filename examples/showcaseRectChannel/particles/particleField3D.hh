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

#ifndef PARTICLE_FIELD_3D_HH
#define PARTICLE_FIELD_3D_HH

#include <utility>

#include "core/globalDefs.h"
#include "particles/particleField3D.h"

namespace plb {

/* *************** class ParticleField3D ************************************ */

template <typename T, template <typename U> class Descriptor>
ParticleField3D<T, Descriptor>::ParticleField3D(
    plint nx, plint ny, plint nz, BlockDataTransfer3D *dataTransfer) :
    AtomicBlock3D(nx, ny, nz, dataTransfer)
{ }

template <typename T, template <typename U> class Descriptor>
bool ParticleField3D<T, Descriptor>::isContained(Array<T, 3> const &particlePos, Box3D box) const
{
    Dot3D const &location = this->getLocation();
    T x = particlePos[0] - location.x;
    T y = particlePos[1] - location.y;
    T z = particlePos[2] - location.z;

    return (x > (T)box.x0 - (T)0.5) && (x <= (T)box.x1 + (T)0.5) && (y > (T)box.y0 - (T)0.5)
           && (y <= (T)box.y1 + (T)0.5) && (z > (T)box.z0 - (T)0.5) && (z <= (T)box.z1 + (T)0.5);
}

template <typename T, template <typename U> class Descriptor>
plint ParticleField3D<T, Descriptor>::nearestCell(T pos)
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
void ParticleField3D<T, Descriptor>::computeGridPosition(
    Array<T, 3> const &position, plint &iX, plint &iY, plint &iZ) const
{
    Dot3D const &location = this->getLocation();
    iX = nearestCell(position[0]) - location.x;
    iY = nearestCell(position[1]) - location.y;
    iZ = nearestCell(position[2]) - location.z;
}

/* *************** class DenseParticleDataTransfer3D ************************ */

template <typename T, template <typename U> class Descriptor>
DenseParticleDataTransfer3D<T, Descriptor>::DenseParticleDataTransfer3D() :
    particleField(0), constParticleField(0)
{ }

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer3D<T, Descriptor>::setBlock(AtomicBlock3D &block)
{
    particleField = dynamic_cast<DenseParticleField3D<T, Descriptor> *>(&block);
    PLB_ASSERT(particleField);
    constParticleField = particleField;
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer3D<T, Descriptor>::setConstBlock(AtomicBlock3D const &block)
{
    constParticleField = dynamic_cast<DenseParticleField3D<T, Descriptor> const *>(&block);
    PLB_ASSERT(constParticleField);
}

template <typename T, template <typename U> class Descriptor>
DenseParticleDataTransfer3D<T, Descriptor> *DenseParticleDataTransfer3D<T, Descriptor>::clone()
    const
{
    return new DenseParticleDataTransfer3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
plint DenseParticleDataTransfer3D<T, Descriptor>::staticCellSize() const
{
    return 0;  // Particle containers have only dynamic data.
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer3D<T, Descriptor>::send(
    Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    PLB_PRECONDITION(constParticleField);
    buffer.clear();
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the send procedure
    //   is run whenever kind is one of the dynamic types.
    if ((kind == modif::dynamicVariables) || (kind == modif::allVariables)
        || (kind == modif::dataStructure))
    {
        std::vector<Particle3D<T, Descriptor> const *> foundParticles;
        constParticleField->findParticles(domain, foundParticles);
        for (pluint iParticle = 0; iParticle < foundParticles.size(); ++iParticle) {
            // The serialize function automatically reallocates memory for buffer.
            serialize(*foundParticles[iParticle], buffer);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer3D<T, Descriptor>::receive(
    Box3D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(particleField);
    PLB_PRECONDITION(contained(domain, particleField->getBoundingBox()));
    // Clear the existing data before introducing the new data.
    particleField->removeParticles(domain);
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
            Particle3D<T, Descriptor> *newParticle =
                meta::particleRegistration3D<T, Descriptor>().generate(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            particleField->addParticle(domain, newParticle);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer3D<T, Descriptor>::receive(
    Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset)
{
    PLB_PRECONDITION(particleField);
    if (absoluteOffset.x == 0 && absoluteOffset.y == 0 && absoluteOffset.z == 0) {
        receive(domain, buffer, kind);
        return;
    }
    PLB_PRECONDITION(contained(domain, particleField->getBoundingBox()));
    Array<T, 3> realAbsoluteOffset((T)absoluteOffset.x, (T)absoluteOffset.y, (T)absoluteOffset.z);
    // Clear the existing data before introducing the new data.
    particleField->removeParticles(domain);
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
            Particle3D<T, Descriptor> *newParticle =
                meta::particleRegistration3D<T, Descriptor>().generate(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            newParticle->getPosition() += realAbsoluteOffset;

            particleField->addParticle(domain, newParticle);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer3D<T, Descriptor>::attribute(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
    modif::ModifT kind)
{
    Box3D fromDomain(toDomain.shift(deltaX, deltaY, deltaZ));
    std::vector<char> buffer;
    DenseParticleField3D<T, Descriptor> const &fromParticleField =
        dynamic_cast<DenseParticleField3D<T, Descriptor> const &>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind);
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleDataTransfer3D<T, Descriptor>::attribute(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
    modif::ModifT kind, Dot3D absoluteOffset)
{
    Box3D fromDomain(toDomain.shift(deltaX, deltaY, deltaZ));
    std::vector<char> buffer;
    DenseParticleField3D<T, Descriptor> const &fromParticleField =
        dynamic_cast<DenseParticleField3D<T, Descriptor> const &>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind, absoluteOffset);
}

/* *************** class DenseParticleField3D ********************** */

template <typename T, template <typename U> class Descriptor>
DenseParticleField3D<T, Descriptor>::DenseParticleField3D(plint nx, plint ny, plint nz) :
    ParticleField3D<T, Descriptor>(nx, ny, nz, new DenseParticleDataTransfer3D<T, Descriptor>()),
    particleGrid(nx, ny, nz)
{ }

template <typename T, template <typename U> class Descriptor>
DenseParticleField3D<T, Descriptor>::~DenseParticleField3D()
{
    for (plint iX = 0; iX < particleGrid.getNx(); ++iX) {
        for (plint iY = 0; iY < particleGrid.getNy(); ++iY) {
            for (plint iZ = 0; iZ < particleGrid.getNz(); ++iZ) {
                for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY, iZ).size();
                     ++iParticle) {
                    delete particleGrid.get(iX, iY, iZ)[iParticle];
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
DenseParticleField3D<T, Descriptor>::DenseParticleField3D(DenseParticleField3D const &rhs) :
    ParticleField3D<T, Descriptor>(rhs),
    particleGrid(rhs.particleGrid.getNx(), rhs.particleGrid.getNy(), rhs.particleGrid.getNz())
{
    for (plint iX = 0; iX < particleGrid.getNx(); ++iX) {
        for (plint iY = 0; iY < particleGrid.getNy(); ++iY) {
            for (plint iZ = 0; iZ < particleGrid.getNz(); ++iZ) {
                particleGrid.get(iX, iY, iZ).resize(rhs.particleGrid.get(iX, iY, iZ).size());
                for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY, iZ).size();
                     ++iParticle) {
                    particleGrid.get(iX, iY, iZ)[iParticle] =
                        rhs.particleGrid.get(iX, iY, iZ)[iParticle]->clone();
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
DenseParticleField3D<T, Descriptor> &DenseParticleField3D<T, Descriptor>::operator=(
    DenseParticleField3D<T, Descriptor> const &rhs)
{
    DenseParticleField3D<T, Descriptor>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
DenseParticleField3D<T, Descriptor> *DenseParticleField3D<T, Descriptor>::clone() const
{
    return new DenseParticleField3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::swap(DenseParticleField3D<T, Descriptor> &rhs)
{
    ParticleField3D<T, Descriptor>::swap(rhs);
    particleGrid.swap(rhs.particleGrid);
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::addParticle(
    Box3D domain, Particle3D<T, Descriptor> *particle)
{
    plint iX, iY, iZ;
    this->computeGridPosition(particle->getPosition(), iX, iY, iZ);
    Box3D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)
        && contained(iX, iY, iZ, finalDomain))
    {
        particleGrid.get(iX, iY, iZ).push_back(particle);
    } else {
        delete particle;
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::removeParticles(Box3D domain)
{
    Box3D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (plint iZ = finalDomain.z0; iZ <= finalDomain.z1; ++iZ) {
                    for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY, iZ).size();
                         ++iParticle) {
                        delete particleGrid.get(iX, iY, iZ)[iParticle];
                    }
                    particleGrid.get(iX, iY, iZ).clear();
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::removeParticles(Box3D domain, plint tag)
{
    Box3D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (plint iZ = finalDomain.z0; iZ <= finalDomain.z1; ++iZ) {
                    typename std::vector<Particle3D<T, Descriptor> *>::iterator it =
                        particleGrid.get(iX, iY, iZ).begin();
                    for (; it != particleGrid.get(iX, iY, iZ).end();) {
                        if ((*it)->getTag() == tag) {
                            delete *it;
                            it = particleGrid.get(iX, iY, iZ).erase(it);
                        } else {
                            ++it;
                        }
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::findParticles(
    Box3D domain, std::vector<Particle3D<T, Descriptor> *> &found)
{
    found.clear();
    PLB_ASSERT(contained(domain, particleGrid.getBoundingBox()));
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY, iZ).size();
                     ++iParticle) {
                    found.push_back(particleGrid.get(iX, iY, iZ)[iParticle]);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::findParticles(
    Box3D domain, std::vector<Particle3D<T, Descriptor> const *> &found) const
{
    found.clear();
    PLB_ASSERT(contained(domain, particleGrid.getBoundingBox()));
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY, iZ).size();
                     ++iParticle) {
                    found.push_back(particleGrid.get(iX, iY, iZ)[iParticle]);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::velocityToParticleCoupling(
    Box3D domain, TensorField3D<T, 3> &velocityField, T scaling)
{
    Box3D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (plint iZ = finalDomain.z0; iZ <= finalDomain.z1; ++iZ) {
                    for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY, iZ).size();
                         ++iParticle) {
                        particleGrid.get(iX, iY, iZ)[iParticle]->velocityToParticle(
                            velocityField, scaling);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::velocityToParticleCoupling(
    Box3D domain, NTensorField3D<T> &velocityField, T scaling)
{
    Box3D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (plint iZ = finalDomain.z0; iZ <= finalDomain.z1; ++iZ) {
                    for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY, iZ).size();
                         ++iParticle) {
                        particleGrid.get(iX, iY, iZ)[iParticle]->velocityToParticle(
                            velocityField, scaling);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::rhoBarJtoParticleCoupling(
    Box3D domain, NTensorField3D<T> &rhoBarJfield, bool velIsJ, T scaling)
{
    Box3D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (plint iZ = finalDomain.z0; iZ <= finalDomain.z1; ++iZ) {
                    for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY, iZ).size();
                         ++iParticle) {
                        particleGrid.get(iX, iY, iZ)[iParticle]->rhoBarJtoParticle(
                            rhoBarJfield, velIsJ, scaling);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::fluidToParticleCoupling(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, T scaling)
{
    Box3D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (plint iZ = finalDomain.z0; iZ <= finalDomain.z1; ++iZ) {
                    for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY, iZ).size();
                         ++iParticle) {
                        particleGrid.get(iX, iY, iZ)[iParticle]->fluidToParticle(lattice, scaling);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::fluidToParticleCoupling(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &sF, T scaling)
{
    Box3D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (plint iZ = finalDomain.z0; iZ <= finalDomain.z1; ++iZ) {
                    for (pluint iParticle = 0; iParticle < particleGrid.get(iX, iY, iZ).size();
                         ++iParticle) {
                        particleGrid.get(iX, iY, iZ)[iParticle]->fluidToParticle(
                            lattice, sF, scaling);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void DenseParticleField3D<T, Descriptor>::advanceParticles(Box3D domain, T cutOffSpeedSqr)
{
    Box3D finalDomain;
    if (intersect(domain, particleGrid.getBoundingBox(), finalDomain)) {
        std::vector<std::pair<Dot3D, Particle3D<T, Descriptor> *> > nextCellParticles;
        for (plint iX = finalDomain.x0; iX <= finalDomain.x1; ++iX) {
            for (plint iY = finalDomain.y0; iY <= finalDomain.y1; ++iY) {
                for (plint iZ = finalDomain.z0; iZ <= finalDomain.z1; ++iZ) {
                    std::vector<Particle3D<T, Descriptor> *> &particles =
                        particleGrid.get(iX, iY, iZ);
                    std::vector<Particle3D<T, Descriptor> *> newLocalParticles;
                    for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
                        Particle3D<T, Descriptor> *particle = particles[iParticle];
                        Array<T, 3> oldPos(particle->getPosition());
                        particle->advance();
                        if (cutOffSpeedSqr >= T()
                            && normSqr(oldPos - particle->getPosition()) < cutOffSpeedSqr) {
                            delete particle;
                        } else {
                            plint newX, newY, newZ;
                            this->computeGridPosition(particle->getPosition(), newX, newY, newZ);
                            if (newX == iX && newY == iY && newZ == iZ) {
                                newLocalParticles.push_back(particle);
                            } else {
                                if (contained(newX, newY, newZ, finalDomain)) {
                                    nextCellParticles.push_back(
                                        std::make_pair(Dot3D(newX, newY, newZ), particle));
                                } else {
                                    delete particle;
                                }
                            }
                        }
                    }
                    newLocalParticles.swap(particles);
                }
            }
        }
        for (pluint i = 0; i < nextCellParticles.size(); ++i) {
            plint newX = nextCellParticles[i].first.x;
            plint newY = nextCellParticles[i].first.y;
            plint newZ = nextCellParticles[i].first.z;
            particleGrid.get(newX, newY, newZ).push_back(nextCellParticles[i].second);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
std::string DenseParticleField3D<T, Descriptor>::getBlockName()
{
    return std::string("DenseParticleField3D");
}

template <typename T, template <typename U> class Descriptor>
std::string DenseParticleField3D<T, Descriptor>::basicType()
{
    return std::string(NativeType<T>::getName());
}

template <typename T, template <typename U> class Descriptor>
std::string DenseParticleField3D<T, Descriptor>::descriptorType()
{
    return std::string(Descriptor<T>::name);
}

/* *************** class LightParticleDataTransfer3D ************************ */

template <typename T, template <typename U> class Descriptor>
LightParticleDataTransfer3D<T, Descriptor>::LightParticleDataTransfer3D() :
    particleField(0), constParticleField(0)
{ }

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer3D<T, Descriptor>::setBlock(AtomicBlock3D &block)
{
    particleField = dynamic_cast<LightParticleField3D<T, Descriptor> *>(&block);
    PLB_ASSERT(particleField);
    constParticleField = particleField;
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer3D<T, Descriptor>::setConstBlock(AtomicBlock3D const &block)
{
    constParticleField = dynamic_cast<LightParticleField3D<T, Descriptor> const *>(&block);
    PLB_ASSERT(constParticleField);
}

template <typename T, template <typename U> class Descriptor>
LightParticleDataTransfer3D<T, Descriptor> *LightParticleDataTransfer3D<T, Descriptor>::clone()
    const
{
    return new LightParticleDataTransfer3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
plint LightParticleDataTransfer3D<T, Descriptor>::staticCellSize() const
{
    return 0;  // Particle containers have only dynamic data.
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer3D<T, Descriptor>::send(
    Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    PLB_PRECONDITION(constParticleField);
    buffer.clear();
    // Particles, by definition, are dynamic data, and they need to
    //   be reconstructed in any case. Therefore, the send procedure
    //   is run whenever kind is one of the dynamic types.
    if ((kind == modif::dynamicVariables) || (kind == modif::allVariables)
        || (kind == modif::dataStructure))
    {
        std::vector<Particle3D<T, Descriptor> const *> foundParticles;
        constParticleField->findParticles(domain, foundParticles);
        for (pluint iParticle = 0; iParticle < foundParticles.size(); ++iParticle) {
            // The serialize function automatically reallocates memory for buffer.
            serialize(*foundParticles[iParticle], buffer);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer3D<T, Descriptor>::receive(
    Box3D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(particleField);
    PLB_PRECONDITION(contained(domain, particleField->getBoundingBox()));
    // Clear the existing data before introducing the new data.
    particleField->removeParticles(domain);
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
            Particle3D<T, Descriptor> *newParticle =
                meta::particleRegistration3D<T, Descriptor>().generate(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            particleField->addParticle(domain, newParticle);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer3D<T, Descriptor>::receive(
    Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset)
{
    PLB_PRECONDITION(particleField);
    if (absoluteOffset.x == 0 && absoluteOffset.y == 0 && absoluteOffset.z == 0) {
        receive(domain, buffer, kind);
        return;
    }
    PLB_PRECONDITION(contained(domain, particleField->getBoundingBox()));
    Array<T, 3> realAbsoluteOffset((T)absoluteOffset.x, (T)absoluteOffset.y, (T)absoluteOffset.z);
    // Clear the existing data before introducing the new data.
    particleField->removeParticles(domain);
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
            Particle3D<T, Descriptor> *newParticle =
                meta::particleRegistration3D<T, Descriptor>().generate(unserializer);
            posInBuffer = unserializer.getCurrentPos();
            newParticle->getPosition() += realAbsoluteOffset;

            particleField->addParticle(domain, newParticle);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer3D<T, Descriptor>::attribute(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
    modif::ModifT kind)
{
    Box3D fromDomain(toDomain.shift(deltaX, deltaY, deltaZ));
    std::vector<char> buffer;
    LightParticleField3D<T, Descriptor> const &fromParticleField =
        dynamic_cast<LightParticleField3D<T, Descriptor> const &>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind);
}

template <typename T, template <typename U> class Descriptor>
void LightParticleDataTransfer3D<T, Descriptor>::attribute(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
    modif::ModifT kind, Dot3D absoluteOffset)
{
    Box3D fromDomain(toDomain.shift(deltaX, deltaY, deltaZ));
    std::vector<char> buffer;
    LightParticleField3D<T, Descriptor> const &fromParticleField =
        dynamic_cast<LightParticleField3D<T, Descriptor> const &>(from);
    fromParticleField.getDataTransfer().send(fromDomain, buffer, kind);
    receive(toDomain, buffer, kind, absoluteOffset);
}

/* *************** class LightParticleField3D ********************** */

template <typename T, template <typename U> class Descriptor>
LightParticleField3D<T, Descriptor>::LightParticleField3D(plint nx, plint ny, plint nz) :
    ParticleField3D<T, Descriptor>(nx, ny, nz, new LightParticleDataTransfer3D<T, Descriptor>())
{ }

template <typename T, template <typename U> class Descriptor>
LightParticleField3D<T, Descriptor>::~LightParticleField3D()
{
    for (pluint i = 0; i < particles.size(); ++i) {
        delete particles[i];
    }
}

template <typename T, template <typename U> class Descriptor>
LightParticleField3D<T, Descriptor>::LightParticleField3D(LightParticleField3D const &rhs) :
    ParticleField3D<T, Descriptor>(rhs)
{
    for (pluint i = 0; i < rhs.particles.size(); ++i) {
        particles.push_back(rhs.particles[i]->clone());
    }
}

template <typename T, template <typename U> class Descriptor>
LightParticleField3D<T, Descriptor> &LightParticleField3D<T, Descriptor>::operator=(
    LightParticleField3D<T, Descriptor> const &rhs)
{
    LightParticleField3D<T, Descriptor>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
LightParticleField3D<T, Descriptor> *LightParticleField3D<T, Descriptor>::clone() const
{
    return new LightParticleField3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField3D<T, Descriptor>::swap(LightParticleField3D<T, Descriptor> &rhs)
{
    ParticleField3D<T, Descriptor>::swap(rhs);
    particles.swap(rhs.particles);
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField3D<T, Descriptor>::addParticle(
    Box3D domain, Particle3D<T, Descriptor> *particle)
{
    Box3D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)
        && this->isContained(particle->getPosition(), finalDomain))
    {
        particles.push_back(particle);
    } else {
        delete particle;
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField3D<T, Descriptor>::removeParticles(Box3D domain)
{
    std::vector<Particle3D<T, Descriptor> *> remainingParticles;
    Box3D finalDomain;
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
void LightParticleField3D<T, Descriptor>::removeParticles(Box3D domain, plint tag)
{
    std::vector<Particle3D<T, Descriptor> *> remainingParticles;
    Box3D finalDomain;
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
void LightParticleField3D<T, Descriptor>::findParticles(
    Box3D domain, std::vector<Particle3D<T, Descriptor> *> &found)
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
void LightParticleField3D<T, Descriptor>::findParticles(
    Box3D domain, std::vector<Particle3D<T, Descriptor> const *> &found) const
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
void LightParticleField3D<T, Descriptor>::velocityToParticleCoupling(
    Box3D domain, TensorField3D<T, 3> &velocityField, T scaling)
{
    Box3D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            if (this->isContained(particles[i]->getPosition(), finalDomain)) {
                particles[i]->velocityToParticle(velocityField, scaling);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField3D<T, Descriptor>::velocityToParticleCoupling(
    Box3D domain, NTensorField3D<T> &velocityField, T scaling)
{
    Box3D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            if (this->isContained(particles[i]->getPosition(), finalDomain)) {
                particles[i]->velocityToParticle(velocityField, scaling);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField3D<T, Descriptor>::rhoBarJtoParticleCoupling(
    Box3D domain, NTensorField3D<T> &rhoBarJfield, bool velIsJ, T scaling)
{
    Box3D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            if (this->isContained(particles[i]->getPosition(), finalDomain)) {
                particles[i]->rhoBarJtoParticle(rhoBarJfield, velIsJ, scaling);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField3D<T, Descriptor>::fluidToParticleCoupling(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, T scaling)
{
    Box3D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            if (this->isContained(particles[i]->getPosition(), finalDomain)) {
                particles[i]->fluidToParticle(lattice, scaling);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void LightParticleField3D<T, Descriptor>::advanceParticles(Box3D domain, T cutOffSpeedSqr)
{
    std::vector<Particle3D<T, Descriptor> *> remainingParticles;
    Box3D finalDomain;
    if (intersect(domain, this->getBoundingBox(), finalDomain)) {
        for (pluint i = 0; i < particles.size(); ++i) {
            Particle3D<T, Descriptor> *particle = particles[i];
            if (this->isContained(particle->getPosition(), finalDomain)) {
                Array<T, 3> oldPos(particle->getPosition());
                particle->advance();
                if ((cutOffSpeedSqr >= T()
                     && normSqr(oldPos - particle->getPosition()) < cutOffSpeedSqr)
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
std::string LightParticleField3D<T, Descriptor>::getBlockName()
{
    return std::string("LightParticleField3D");
}

template <typename T, template <typename U> class Descriptor>
std::string LightParticleField3D<T, Descriptor>::basicType()
{
    return std::string(NativeType<T>::getName());
}

template <typename T, template <typename U> class Descriptor>
std::string LightParticleField3D<T, Descriptor>::descriptorType()
{
    return std::string(Descriptor<T>::name);
}

}  // namespace plb

#endif  // PARTICLE_FIELD_3D_HH
