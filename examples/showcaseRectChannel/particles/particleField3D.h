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

#ifndef PARTICLE_FIELD_3D_H
#define PARTICLE_FIELD_3D_H

#include <vector>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/globalDefs.h"
#include "particles/particle3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class ParticleField3D : public AtomicBlock3D {
public:
    ParticleField3D(plint nx, plint ny, plint nz, BlockDataTransfer3D *dataTransfer);
    /// Add a particle if it is part of the domain, else delete it.
    /** This method with domain-argument is provided here exclusively,
     *    because it may not be easy for an outside instance to decide
     *    autonomously whether a particle is inside a domain or not
     *    (because domains are enlarged by half a lattice site).
     **/
    virtual void addParticle(Box3D domain, Particle3D<T, Descriptor> *particle) = 0;
    /// Remove all particles found in the indicated domain.
    virtual void removeParticles(Box3D domain) = 0;
    /// Remove all particles of a certain tag found in the indicated domain.
    virtual void removeParticles(Box3D domain, plint tag) = 0;
    /// Return all particles found in the indicated domain.
    virtual void findParticles(Box3D domain, std::vector<Particle3D<T, Descriptor> *> &found) = 0;
    /// Return all particles found in the indicated domain; const version
    virtual void findParticles(
        Box3D domain, std::vector<Particle3D<T, Descriptor> const *> &found) const = 0;
    /// Execute velocity-particle interaction for all particles contained in the domain.
    virtual void velocityToParticleCoupling(
        Box3D domain, TensorField3D<T, 3> &velocity, T scaling = 0.) = 0;
    /// Execute velocity-particle interaction for all particles contained in the domain.
    virtual void velocityToParticleCoupling(
        Box3D domain, NTensorField3D<T> &velocity, T scaling = 0.) = 0;
    /// Execute velocity-particle interaction for all particles contained in the domain.
    virtual void rhoBarJtoParticleCoupling(
        Box3D domain, NTensorField3D<T> &rhoBarJ, bool velIsJ, T scaling = 0.) = 0;
    /// Execute fluid-particle interaction for all particles contained in the domain.
    virtual void fluidToParticleCoupling(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, T scaling = 0.) = 0;
    /// Advance all particles contained in the domain. When the speed of a particle drops
    ///   below sqrt(cutOffSpeedSqr), the particle is eliminated. Negative cutOffSpeedSqr means
    ///   no cutoff.
    virtual void advanceParticles(Box3D domain, T cutOffSpeedSqr = -1.) = 0;
    virtual identifiers::BlockId getBlockId() const
    {
        return identifiers::ParticleId;
    }

public:
    /// Helper function: returns if a given particle is situated in the indicated (local) domain.
    bool isContained(Array<T, 3> const &particlePos, Box3D box) const;
    void computeGridPosition(Array<T, 3> const &position, plint &iX, plint &iY, plint &iZ) const;

private:
    static plint nearestCell(T pos);
};

template <typename T, template <typename U> class Descriptor>
class DenseParticleField3D;

template <typename T, template <typename U> class Descriptor>
class DenseParticleDataTransfer3D : public BlockDataTransfer3D {
public:
    DenseParticleDataTransfer3D();
    virtual void setBlock(AtomicBlock3D &block);
    virtual void setConstBlock(AtomicBlock3D const &block);
    virtual DenseParticleDataTransfer3D<T, Descriptor> *clone() const;
    virtual plint staticCellSize() const;
    virtual void send(Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    virtual void receive(Box3D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset);
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    {
        receive(domain, buffer, kind);
    }
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind);
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind, Dot3D absoluteOffset);

private:
    DenseParticleField3D<T, Descriptor> *particleField;
    DenseParticleField3D<T, Descriptor> const *constParticleField;
};

template <typename T, template <typename U> class Descriptor>
class DenseParticleField3D : public ParticleField3D<T, Descriptor> {
public:
    typedef Particle3D<T, Descriptor> ParticleT;

public:
    DenseParticleField3D(plint nx, plint ny, plint nz);
    virtual ~DenseParticleField3D();
    DenseParticleField3D(DenseParticleField3D<T, Descriptor> const &rhs);
    DenseParticleField3D<T, Descriptor> &operator=(DenseParticleField3D<T, Descriptor> const &rhs);
    DenseParticleField3D<T, Descriptor> *clone() const;
    void swap(DenseParticleField3D<T, Descriptor> &rhs);

public:
    virtual void addParticle(Box3D domain, Particle3D<T, Descriptor> *particle);
    virtual void removeParticles(Box3D domain);
    virtual void removeParticles(Box3D domain, plint tag);
    virtual void findParticles(Box3D domain, std::vector<Particle3D<T, Descriptor> *> &found);
    virtual void findParticles(
        Box3D domain, std::vector<Particle3D<T, Descriptor> const *> &found) const;
    virtual void velocityToParticleCoupling(
        Box3D domain, TensorField3D<T, 3> &velocity, T scaling = 0.);
    virtual void velocityToParticleCoupling(
        Box3D domain, NTensorField3D<T> &velocity, T scaling = 0.);
    virtual void rhoBarJtoParticleCoupling(
        Box3D domain, NTensorField3D<T> &rhoBarJ, bool velIsJ, T scaling = 0.);
    virtual void fluidToParticleCoupling(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, T scaling = 0.);
    virtual void fluidToParticleCoupling(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &sF, T scaling = 0.);
    virtual void advanceParticles(Box3D domain, T cutOffSpeedSqr = -1.);

public:
    static std::string getBlockName();
    static std::string basicType();
    static std::string descriptorType();

private:
    ScalarField3D<std::vector<Particle3D<T, Descriptor> *> > particleGrid;
    DenseParticleDataTransfer3D<T, Descriptor> dataTransfer;
};

template <typename T, template <typename U> class Descriptor>
class LightParticleField3D;

template <typename T, template <typename U> class Descriptor>
class LightParticleDataTransfer3D : public BlockDataTransfer3D {
public:
    LightParticleDataTransfer3D();
    virtual void setBlock(AtomicBlock3D &block);
    virtual void setConstBlock(AtomicBlock3D const &block);
    virtual LightParticleDataTransfer3D<T, Descriptor> *clone() const;
    virtual plint staticCellSize() const;
    virtual void send(Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    virtual void receive(Box3D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset);
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    {
        receive(domain, buffer, kind);
    }
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind);
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind, Dot3D absoluteOffset);

private:
    LightParticleField3D<T, Descriptor> *particleField;
    LightParticleField3D<T, Descriptor> const *constParticleField;
};

template <typename T, template <typename U> class Descriptor>
class LightParticleField3D : public ParticleField3D<T, Descriptor> {
public:
    typedef Particle3D<T, Descriptor> ParticleT;

public:
    LightParticleField3D(plint nx, plint ny, plint nz);
    virtual ~LightParticleField3D();
    LightParticleField3D(LightParticleField3D<T, Descriptor> const &rhs);
    LightParticleField3D<T, Descriptor> &operator=(LightParticleField3D<T, Descriptor> const &rhs);
    LightParticleField3D<T, Descriptor> *clone() const;
    void swap(LightParticleField3D<T, Descriptor> &rhs);

public:
    virtual void addParticle(Box3D domain, Particle3D<T, Descriptor> *particle);
    virtual void removeParticles(Box3D domain);
    virtual void removeParticles(Box3D domain, plint tag);
    virtual void findParticles(Box3D domain, std::vector<Particle3D<T, Descriptor> *> &found);
    virtual void findParticles(
        Box3D domain, std::vector<Particle3D<T, Descriptor> const *> &found) const;
    virtual void velocityToParticleCoupling(
        Box3D domain, TensorField3D<T, 3> &velocity, T scaling = 0.);
    virtual void velocityToParticleCoupling(
        Box3D domain, NTensorField3D<T> &velocity, T scaling = 0.);
    virtual void rhoBarJtoParticleCoupling(
        Box3D domain, NTensorField3D<T> &rhoBarJ, bool velIsJ, T scaling = 0.);
    virtual void fluidToParticleCoupling(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, T scaling = 0.);
    virtual void advanceParticles(Box3D domain, T cutOffSpeedSqr = -1.);

public:
    static std::string getBlockName();
    static std::string basicType();
    static std::string descriptorType();

private:
    std::vector<Particle3D<T, Descriptor> *> particles;
    LightParticleDataTransfer3D<T, Descriptor> dataTransfer;
};

}  // namespace plb

#endif  // PARTICLE_FIELD_3D_H
