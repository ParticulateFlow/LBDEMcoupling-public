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

#ifndef PARTICLE_FIELD_2D_H
#define PARTICLE_FIELD_2D_H

#include <vector>

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/dataField2D.h"
#include "core/globalDefs.h"
#include "particles/particle2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class ParticleField2D : public AtomicBlock2D {
public:
    ParticleField2D(plint nx, plint ny);
    /// Add a particle if it is part of the domain, else delete it.
    /** This method with domain-argument is provided here exclusively,
     *    because it may not be easy for an outside instance to decide
     *    autonomously whether a particle is inside a domain or not
     *    (because domains are enlarged by half a lattice site).
     **/
    virtual void addParticle(Box2D domain, Particle2D<T, Descriptor> *particle) = 0;
    /// Remove all particles found in the indicated domain.
    virtual void removeParticles(Box2D domain) = 0;
    /// Remove all particles of a certain tag found in the indicated domain.
    virtual void removeParticles(Box2D domain, plint tag) = 0;
    /// Return all particles found in the indicated domain.
    virtual void findParticles(Box2D domain, std::vector<Particle2D<T, Descriptor> *> &found) = 0;
    /// Return all particles found in the indicated domain; const version
    virtual void findParticles(
        Box2D domain, std::vector<Particle2D<T, Descriptor> const *> &found) const = 0;
    /// Execute velocity-particle interaction for all particles contained in the domain.
    virtual void velocityToParticleCoupling(
        Box2D domain, TensorField2D<T, 2> &velocity, T scaling = 0.) = 0;
    /// Execute velocity-particle interaction for all particles contained in the domain.
    virtual void velocityToParticleCoupling(
        Box2D domain, NTensorField2D<T> &velocity, T scaling = 0.) = 0;
    /// Execute velocity-particle interaction for all particles contained in the domain.
    virtual void rhoBarJtoParticleCoupling(
        Box2D domain, NTensorField2D<T> &rhoBarJ, bool velIsJ, T scaling = 0.) = 0;
    /// Execute fluid-particle interaction for all particles contained in the domain.
    virtual void fluidToParticleCoupling(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, T scaling = 0.) = 0;
    /// Advance all particles contained in the domain. When the speed of a particle drops
    ///   below sqrt(cutOffValue), the particle is eliminated. Negative cutOffValue means
    ///   no cutoff.
    virtual void advanceParticles(Box2D domain, T cutOffValue = -1.) = 0;
    virtual identifiers::BlockId getBlockId() const
    {
        return identifiers::ParticleId;
    }

public:
    /// Helper function: returns if a given particle is situated in the indicated domain.
    bool isContained(Array<T, 2> const &particlePos, Box2D box) const;
    void computeGridPosition(Array<T, 2> const &position, plint &iX, plint &iY) const;

private:
    static plint nearestCell(T pos);
};

template <typename T, template <typename U> class Descriptor>
class DenseParticleField2D;

template <typename T, template <typename U> class Descriptor>
class DenseParticleDataTransfer2D : public BlockDataTransfer2D {
public:
    DenseParticleDataTransfer2D(DenseParticleField2D<T, Descriptor> &particleField_);
    virtual plint staticCellSize() const;
    virtual void send(Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    virtual void receive(Box2D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot2D absoluteOffset);
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    {
        receive(domain, buffer, kind);
    }
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind);
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind,
        Dot2D absoluteOffset);

private:
    DenseParticleField2D<T, Descriptor> &particleField;
};

template <typename T, template <typename U> class Descriptor>
class DenseParticleField2D : public ParticleField2D<T, Descriptor> {
public:
    typedef Particle2D<T, Descriptor> ParticleT;

public:
    DenseParticleField2D(plint nx, plint ny);
    virtual ~DenseParticleField2D();
    DenseParticleField2D(DenseParticleField2D<T, Descriptor> const &rhs);
    DenseParticleField2D<T, Descriptor> &operator=(DenseParticleField2D<T, Descriptor> const &rhs);
    DenseParticleField2D<T, Descriptor> *clone() const;
    void swap(DenseParticleField2D<T, Descriptor> &rhs);

public:
    virtual void addParticle(Box2D domain, Particle2D<T, Descriptor> *particle);
    virtual void removeParticles(Box2D domain);
    virtual void removeParticles(Box2D domain, plint tag);
    virtual void findParticles(Box2D domain, std::vector<Particle2D<T, Descriptor> *> &found);
    virtual void findParticles(
        Box2D domain, std::vector<Particle2D<T, Descriptor> const *> &found) const;
    virtual void velocityToParticleCoupling(
        Box2D domain, TensorField2D<T, 2> &velocity, T scaling = 0.);
    virtual void velocityToParticleCoupling(
        Box2D domain, NTensorField2D<T> &velocity, T scaling = 0.);
    virtual void rhoBarJtoParticleCoupling(
        Box2D domain, NTensorField2D<T> &rhoBarJ, bool velIsJ, T scaling = 0.);
    virtual void fluidToParticleCoupling(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, T scaling);
    virtual void advanceParticles(Box2D domain, T cutOffValue = -1.);

public:
    virtual DenseParticleDataTransfer2D<T, Descriptor> &getDataTransfer();
    virtual DenseParticleDataTransfer2D<T, Descriptor> const &getDataTransfer() const;
    static std::string getBlockName();
    static std::string basicType();
    static std::string descriptorType();

private:
    ScalarField2D<std::vector<Particle2D<T, Descriptor> *> > particleGrid;
    DenseParticleDataTransfer2D<T, Descriptor> dataTransfer;
};

template <typename T, template <typename U> class Descriptor>
class LightParticleField2D;

template <typename T, template <typename U> class Descriptor>
class LightParticleDataTransfer2D : public BlockDataTransfer2D {
public:
    LightParticleDataTransfer2D(LightParticleField2D<T, Descriptor> &particleField_);
    virtual plint staticCellSize() const;
    virtual void send(Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    virtual void receive(Box2D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot2D absoluteOffset);
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    {
        receive(domain, buffer, kind);
    }
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind);
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind,
        Dot2D absoluteOffset);

private:
    LightParticleField2D<T, Descriptor> &particleField;
};

template <typename T, template <typename U> class Descriptor>
class LightParticleField2D : public ParticleField2D<T, Descriptor> {
public:
    typedef Particle2D<T, Descriptor> ParticleT;

public:
    LightParticleField2D(plint nx, plint ny);
    virtual ~LightParticleField2D();
    LightParticleField2D(LightParticleField2D<T, Descriptor> const &rhs);
    LightParticleField2D<T, Descriptor> &operator=(LightParticleField2D<T, Descriptor> const &rhs);
    LightParticleField2D<T, Descriptor> *clone() const;
    void swap(LightParticleField2D<T, Descriptor> &rhs);

public:
    virtual void addParticle(Box2D domain, Particle2D<T, Descriptor> *particle);
    virtual void removeParticles(Box2D domain);
    virtual void removeParticles(Box2D domain, plint tag);
    virtual void findParticles(Box2D domain, std::vector<Particle2D<T, Descriptor> *> &found);
    virtual void findParticles(
        Box2D domain, std::vector<Particle2D<T, Descriptor> const *> &found) const;
    virtual void velocityToParticleCoupling(
        Box2D domain, TensorField2D<T, 2> &velocity, T scaling = 0.);
    virtual void velocityToParticleCoupling(
        Box2D domain, NTensorField2D<T> &velocity, T scaling = 0.);
    virtual void rhoBarJtoParticleCoupling(
        Box2D domain, NTensorField2D<T> &rhoBarJ, bool velIsJ, T scaling = 0.);
    virtual void fluidToParticleCoupling(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, T scaling = 0.);
    virtual void advanceParticles(Box2D domain, T cutOffValue = -1.);

public:
    virtual LightParticleDataTransfer2D<T, Descriptor> &getDataTransfer();
    virtual LightParticleDataTransfer2D<T, Descriptor> const &getDataTransfer() const;
    static std::string getBlockName();
    static std::string basicType();
    static std::string descriptorType();

private:
    std::vector<Particle2D<T, Descriptor> *> particles;
    LightParticleDataTransfer2D<T, Descriptor> dataTransfer;
};

}  // namespace plb

#endif  // PARTICLE_FIELD_2D_H
