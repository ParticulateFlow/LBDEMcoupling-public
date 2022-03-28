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

#ifndef PARTICLE_3D_H
#define PARTICLE_3D_H

#include <vector>

#include "atomicBlock/blockLattice3D.h"
#include "core/array.h"
#include "core/globalDefs.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class Particle3D {
public:
    Particle3D();
    /// The tag does not need to be unique and exists only for the
    ///   user's convenience.
    Particle3D(plint tag_, Array<T, 3> const &position_);
    virtual ~Particle3D() { }
    /// Couple velocity to particle. Requirement: fluid velocity<0.25 in lattice units.
    /// If this condition is violated, the velocity is trimmed to a 0.25 velocity norm.
    virtual void velocityToParticle(TensorField3D<T, 3> &velocityField, T scaling = 1.) = 0;
    virtual void velocityToParticle(NTensorField3D<T> &velocityField, T scaling = 1.) = 0;
    /// Couple velocity to particle. Requirement: fluid velocity<0.25 in lattice units.
    /// If this condition is violated, the velocity is trimmed to a 0.25 velocity norm.
    virtual void rhoBarJtoParticle(
        NTensorField3D<T> &rhoBarJfield, bool velIsJ, T scaling = 1.) = 0;
    /// Couple fluid velocity to particle. Requirement: fluid velocity<0.25 in lattice units.
    /// If this condition is violated, the velocity is trimmed to a 0.25 velocity norm.
    virtual void fluidToParticle(BlockLattice3D<T, Descriptor> &fluid, T scaling = 1.) = 0;
    virtual void fluidToParticle(
        BlockLattice3D<T, Descriptor> &fluid, ScalarField3D<T> &sF, T scaling = 1.) = 0;
    virtual void advance() = 0;
    Array<T, 3> const &getPosition() const
    {
        return position;
    }
    Array<T, 3> &getPosition()
    {
        return position;
    }
    virtual void reset(Array<T, 3> const &position);
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Return a unique ID for this class.
    virtual int getId() const = 0;
    /// Return this particle's tag (not necessarily unique).
    plint getTag() const;
    /// Set a tag for this particle (not necessarily unique).
    void setTag(plint tag_);
    virtual Particle3D<T, Descriptor> *clone() const = 0;
    virtual bool getScalar(plint whichScalar, T &scalar) const;
    virtual bool setScalar(plint whichScalar, T scalar);
    virtual bool getVector(plint whichVector, Array<T, 3> &vector) const;
    virtual bool setVector(plint whichVector, Array<T, 3> const &vector);
    virtual bool getTensor(plint whichVector, Array<T, SymmetricTensorImpl<T, 3>::n> &tensor) const;
    virtual bool setTensor(plint whichVector, Array<T, SymmetricTensorImpl<T, 3>::n> const &tensor);
    virtual bool setScalars(std::vector<T> const &scalars);
    virtual bool setVectors(std::vector<Array<T, 3> > const &vectors);
    virtual bool setTensors(std::vector<Array<T, SymmetricTensorImpl<T, 3>::n> > const &tensors);
    virtual void rescale(int dxScale, int dtScale);

private:
    plint tag;
    Array<T, 3> position;
};

/// Serialize the whole particle chain into a byte-stream, and append it to
///   the vector data.
template <typename T, template <typename U> class Descriptor>
void serialize(Particle3D<T, Descriptor> const &particle, std::vector<char> &data);

/// Unserialize all data into newly generated particle objects.
template <typename T, template <typename U> class Descriptor>
void generateAndUnserializeParticles(
    std::vector<char> const &data, std::vector<Particle3D<T, Descriptor> *> &particle);

template <typename T, template <typename U> class Descriptor>
class PointParticle3D : public Particle3D<T, Descriptor> {
public:
    PointParticle3D();
    PointParticle3D(plint tag_, Array<T, 3> const &position_, Array<T, 3> const &velocity_);
    virtual void velocityToParticle(TensorField3D<T, 3> &velocityField, T scaling = 1.);
    virtual void velocityToParticle(NTensorField3D<T> &velocityField, T scaling = 1.);
    virtual void rhoBarJtoParticle(NTensorField3D<T> &rhoBarJfield, bool velIsJ, T scaling = 1.);
    virtual void fluidToParticle(BlockLattice3D<T, Descriptor> &fluid, T scaling = 1.);
    virtual void fluidToParticle(
        BlockLattice3D<T, Descriptor> &fluid, ScalarField3D<T> &sF, T scaling = 1.);
    virtual void advance();
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual int getId() const;
    virtual void reset(Array<T, 3> const &position);
    virtual PointParticle3D<T, Descriptor> *clone() const;
    virtual bool getVector(plint whichVector, Array<T, 3> &vector) const;
    virtual bool setVectors(std::vector<Array<T, 3> > const &vectors);
    virtual void rescale(int dxScale, int dtScale);
    Array<T, 3> const &getVelocity() const
    {
        return velocity;
    }
    Array<T, 3> &getVelocity()
    {
        return velocity;
    }

private:
    Array<T, 3> velocity;
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class NormedVelocityParticle3D : public PointParticle3D<T, Descriptor> {
public:
    NormedVelocityParticle3D();
    NormedVelocityParticle3D(
        plint tag_, Array<T, 3> const &position_, Array<T, 3> const &velocity_, T fluidUmax_,
        T particleUmax_, T exponent_);
    virtual void advance();
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual NormedVelocityParticle3D<T, Descriptor> *clone() const;
    virtual void rescale(int dxScale, int dtScale);
    virtual int getId() const;

private:
    T fluidUmax, particleUmax, exponent;
    static int id;
};

/* *************** class RestParticle3D ************************************ */

template <typename T, template <typename U> class Descriptor>
class RestParticle3D : public Particle3D<T, Descriptor> {
public:
    RestParticle3D();
    RestParticle3D(plint tag_, Array<T, 3> const &position);
    virtual void velocityToParticle(TensorField3D<T, 3> &velocityField, T scaling = 1.);
    virtual void velocityToParticle(NTensorField3D<T> &velocityField, T scaling = 1.);
    virtual void rhoBarJtoParticle(NTensorField3D<T> &rhoBarJfield, bool velIsJ, T scaling = 1.);
    virtual void fluidToParticle(BlockLattice3D<T, Descriptor> &fluid, T scaling = 1.);
    virtual void fluidToParticle(
        BlockLattice3D<T, Descriptor> &fluid, ScalarField3D<T> &sF, T scaling = 1.);
    virtual void advance();
    virtual int getId() const;
    virtual RestParticle3D<T, Descriptor> *clone() const;
    virtual bool getVector(plint whichVector, Array<T, 3> &vector) const;

private:
    static int id;
};

/* *************** class VerletParticle3D ************************************ */

template <typename T, template <typename U> class Descriptor>
class VerletParticle3D : public Particle3D<T, Descriptor> {
public:
    VerletParticle3D();
    VerletParticle3D(plint tag_, Array<T, 3> const &position);
    virtual void velocityToParticle(TensorField3D<T, 3> &velocityField, T scaling = 1.);
    virtual void velocityToParticle(NTensorField3D<T> &velocityField, T scaling = 1.);
    virtual void rhoBarJtoParticle(NTensorField3D<T> &rhoBarJfield, bool velIsJ, T scaling = 1.);
    virtual void fluidToParticle(BlockLattice3D<T, Descriptor> &fluid, T scaling = 1.);
    virtual void fluidToParticle(
        BlockLattice3D<T, Descriptor> &fluid, ScalarField3D<T> &sF, T scaling = 1.);
    /// Implements "steps 1 and 2" of the Rest algorithm: given
    ///   x(t), v(t), and a(t), it computes v(t+1/2) and x(t+1).
    virtual void advance();
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual int getId() const;
    virtual void reset(Array<T, 3> const &position);
    virtual VerletParticle3D<T, Descriptor> *clone() const;
    virtual bool getVector(plint whichVector, Array<T, 3> &vector) const;
    virtual bool getScalar(plint whichScalar, T &scalar) const;

    Array<T, 3> const &get_v() const
    {
        return v;
    }
    Array<T, 3> const &get_vHalfTime() const
    {
        return vHalfTime;
    }
    Array<T, 3> const &get_a() const
    {
        return a;
    }
    T getFluidCompliance() const
    {
        return fluidCompliance;
    }
    T get_rho() const
    {
        return rho;
    }
    T get_invRho() const
    {
        return invRho;
    }

    void set_v(Array<T, 3> const &v_)
    {
        v = v_;
    }
    void set_vHalfTime(Array<T, 3> const &vHalfTime_)
    {
        vHalfTime = vHalfTime_;
    }
    void set_a(Array<T, 3> const &a_)
    {
        a = a_;
    }
    void setFluidCompliance(T fluidCompliance_)
    {
        fluidCompliance = fluidCompliance_;
    }
    void set_rho(T rho_)
    {
        rho = rho_;
        invRho = 1. / rho;
    }

private:
    Array<T, 3> v, vHalfTime, a;
    T fluidCompliance;
    T rho, invRho;
    static int id;
};

}  // namespace plb

#endif  // PARTICLE_3D_H
