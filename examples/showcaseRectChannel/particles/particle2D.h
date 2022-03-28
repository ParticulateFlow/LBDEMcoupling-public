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

#ifndef PARTICLE_2D_H
#define PARTICLE_2D_H

#include <vector>

#include "atomicBlock/blockLattice2D.h"
#include "core/array.h"
#include "core/globalDefs.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class Particle2D {
public:
    Particle2D();
    /// The tag does not need to be unique and exists only for the
    ///   user's convenience.
    Particle2D(plint tag_, Array<T, 2> const &position_);
    virtual ~Particle2D() { }
    virtual void velocityToParticle(TensorField2D<T, 2> &velocityField, T scaling = 1.) = 0;
    virtual void velocityToParticle(NTensorField2D<T> &velocityField, T scaling = 1.) = 0;
    virtual void rhoBarJtoParticle(
        NTensorField2D<T> &rhoBarJfield, bool velIsJ, T scaling = 1.) = 0;
    virtual void fluidToParticle(BlockLattice2D<T, Descriptor> &fluid, T scaling = 1.) = 0;
    virtual void advance() = 0;
    Array<T, 2> const &getPosition() const
    {
        return position;
    }
    Array<T, 2> &getPosition()
    {
        return position;
    }
    virtual void reset(Array<T, 2> const &position);
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Return a unique ID for this class.
    virtual int getId() const = 0;
    /// Return this particle's tag (not necessarily unique).
    plint getTag() const;
    /// Set a tag for this particle (not necessarily unique).
    void setTag(plint tag_);
    virtual Particle2D<T, Descriptor> *clone() const = 0;
    virtual bool getScalar(plint whichScalar, T &scalar) const;
    virtual bool setScalar(plint whichScalar, T scalar);
    virtual bool getVector(plint whichVector, Array<T, 2> &vector) const;
    virtual bool setVector(plint whichVector, Array<T, 2> const &vector);
    virtual bool getTensor(plint whichVector, Array<T, SymmetricTensorImpl<T, 2>::n> &tensor) const;
    virtual bool setTensor(plint whichVector, Array<T, SymmetricTensorImpl<T, 2>::n> const &tensor);
    virtual bool setScalars(std::vector<T> const &scalars);
    virtual bool setVectors(std::vector<Array<T, 2> > const &vectors);
    virtual bool setTensors(std::vector<Array<T, SymmetricTensorImpl<T, 2>::n> > const &tensors);
    virtual void rescale(int dxScale, int dtScale);

private:
    plint tag;
    Array<T, 2> position;
};

/// Serialize the whole particle chain into a byte-stream, and append it to
///   the vector data.
template <typename T, template <typename U> class Descriptor>
void serialize(Particle2D<T, Descriptor> const &particle, std::vector<char> &data);

/// Unserialize all data into newly generated particle objects.
template <typename T, template <typename U> class Descriptor>
void generateAndUnserializeParticles(
    std::vector<char> const &data, std::vector<Particle2D<T, Descriptor> *> &particle);

template <typename T, template <typename U> class Descriptor>
class PointParticle2D : public Particle2D<T, Descriptor> {
public:
    PointParticle2D();
    PointParticle2D(plint tag_, Array<T, 2> const &position_, Array<T, 2> const &velocity_);
    virtual void velocityToParticle(TensorField2D<T, 2> &velocityField, T scaling = 1.);
    virtual void velocityToParticle(NTensorField2D<T> &velocityField, T scaling = 1.);
    virtual void rhoBarJtoParticle(NTensorField2D<T> &rhoBarJfield, bool velIsJ, T scaling = 1.);
    virtual void fluidToParticle(BlockLattice2D<T, Descriptor> &fluid, T scaling = 1.);
    virtual void advance();
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual int getId() const;
    virtual void reset(Array<T, 2> const &position);
    virtual PointParticle2D<T, Descriptor> *clone() const;
    virtual bool getVector(plint whichVector, Array<T, 2> &vector) const;
    virtual bool setVectors(std::vector<Array<T, 2> > const &vectors);
    virtual void rescale(int dxScale, int dtScale);
    Array<T, 2> const &getVelocity() const
    {
        return velocity;
    }
    Array<T, 2> &getVelocity()
    {
        return velocity;
    }

private:
    Array<T, 2> velocity;
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class NormedVelocityParticle2D : public PointParticle2D<T, Descriptor> {
public:
    NormedVelocityParticle2D();
    NormedVelocityParticle2D(
        plint tag_, Array<T, 2> const &position_, Array<T, 2> const &velocity_, T fluidUmax_,
        T particleUmax_, T exponent_);
    virtual void advance();
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual NormedVelocityParticle2D<T, Descriptor> *clone() const;
    virtual void rescale(int dxScale, int dtScale);
    virtual int getId() const;

private:
    T fluidUmax, particleUmax, exponent;
    static int id;
};

/* *************** class RestParticle2D ************************************ */

template <typename T, template <typename U> class Descriptor>
class RestParticle2D : public Particle2D<T, Descriptor> {
public:
    RestParticle2D();
    RestParticle2D(plint tag_, Array<T, 2> const &position);
    virtual void velocityToParticle(TensorField2D<T, 2> &velocityField, T scaling = 1.);
    virtual void velocityToParticle(NTensorField2D<T> &velocityField, T scaling = 1.);
    virtual void rhoBarJtoParticle(NTensorField2D<T> &rhoBarJfield, bool velIsJ, T scaling = 1.);
    virtual void fluidToParticle(BlockLattice2D<T, Descriptor> &fluid, T scaling = 1.);
    virtual void advance();
    virtual int getId() const;
    virtual RestParticle2D<T, Descriptor> *clone() const;
    virtual bool getVector(plint whichVector, Array<T, 2> &vector) const;

private:
    static int id;
};

/* *************** class VerletParticle2D ************************************ */

template <typename T, template <typename U> class Descriptor>
class VerletParticle2D : public Particle2D<T, Descriptor> {
public:
    VerletParticle2D();
    VerletParticle2D(plint tag_, Array<T, 2> const &position);
    virtual void velocityToParticle(TensorField2D<T, 2> &velocityField, T scaling = 1.);
    virtual void velocityToParticle(NTensorField2D<T> &velocityField, T scaling = 1.);
    virtual void rhoBarJtoParticle(NTensorField2D<T> &rhoBarJfield, bool velIsJ, T scaling = 1.);
    virtual void fluidToParticle(BlockLattice2D<T, Descriptor> &fluid, T scaling = 1.);
    /// Implements "steps 1 and 2" of the Rest algorithm: given
    ///   x(t), v(t), and a(t), it computes v(t+1/2) and x(t+1).
    virtual void advance();
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual int getId() const;
    virtual void reset(Array<T, 2> const &position);
    virtual VerletParticle2D<T, Descriptor> *clone() const;
    virtual bool getVector(plint whichVector, Array<T, 2> &vector) const;
    virtual bool getScalar(plint whichScalar, T &scalar) const;

    Array<T, 2> const &get_v() const
    {
        return v;
    }
    Array<T, 2> const &get_vHalfTime() const
    {
        return vHalfTime;
    }
    Array<T, 2> const &get_a() const
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

    void set_v(Array<T, 2> const &v_)
    {
        v = v_;
    }
    void set_vHalfTime(Array<T, 2> const &vHalfTime_)
    {
        vHalfTime = vHalfTime_;
    }
    void set_a(Array<T, 2> const &a_)
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
    Array<T, 2> v, vHalfTime, a;
    T fluidCompliance;
    T rho, invRho;
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class VerletGeryaParticle2D : public VerletParticle2D<T, Descriptor> {
public:
    VerletGeryaParticle2D();
    VerletGeryaParticle2D(plint tag_, Array<T, 2> const &position_, T fluidConcentration_);
    // virtual void advance();
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual VerletGeryaParticle2D<T, Descriptor> *clone() const;
    virtual int getId() const;
    virtual void reset(Array<T, 2> const &position);

    T getFluidConcentration();

    void setFluidConcentration(T fluidConcentration_);

private:
    T fluidConcentration;
    static int id;
};

}  // namespace plb

#endif  // PARTICLE_2D_H
