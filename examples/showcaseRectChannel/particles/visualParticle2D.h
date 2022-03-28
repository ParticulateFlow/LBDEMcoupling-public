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

#ifndef VISUAL_PARTICLE_2D_H
#define VISUAL_PARTICLE_2D_H

#include <vector>

#include "atomicBlock/blockLattice2D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "particles/particle2D.h"
#include "particles/particleIdentifiers2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class VisualParticle2D : public Particle2D<T, Descriptor> {
public:
    VisualParticle2D();
    VisualParticle2D(plint tag_, Array<T, 2> const &position);
    VisualParticle2D(
        plint tag_, Array<T, 2> const &position, std::vector<T> const &scalars_,
        std::vector<Array<T, 2> > const &vectors_);
    virtual void velocityToParticle(TensorField2D<T, 2> &velocityField, T scaling = 1.) { }
    virtual void velocityToParticle(NTensorField2D<T> &velocityField, T scaling = 1.) { }
    virtual void rhoBarJtoParticle(NTensorField2D<T> &rhoBarJfield, bool velIsJ, T scaling = 1.) { }
    virtual void fluidToParticle(BlockLattice2D<T, Descriptor> &fluid, T scaling = 1.) { }
    virtual void advance() { }
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual int getId() const;
    virtual void reset(Array<T, 2> const &position);
    virtual VisualParticle2D<T, Descriptor> *clone() const;
    virtual bool getVector(plint whichVector, Array<T, 2> &vector) const;
    virtual bool getScalar(plint whichScalar, T &scalar) const;
    virtual bool setScalar(plint whichScalar, T scalar);
    virtual bool setScalars(std::vector<T> const &scalars_);
    virtual bool setVector(plint whichVector, Array<T, 2> const &vector);
    virtual bool setVectors(std::vector<Array<T, 2> > const &vectors_);

private:
    std::vector<T> scalars;
    std::vector<Array<T, 2> > vectors;
    static int id;
};

// template<typename T, template<typename U> class Descriptor>
// class MappingParticle2D : public PointParticle2D<T,Descriptor> {
// public:
//     MappingParticle2D();
//     MappingParticle2D(plint tag_, Array<T,2> const& position_, Array<T,2> const& velocity_);
//     virtual void advance();
//     virtual void serialize(HierarchicSerializer& serializer) const;
//     virtual void unserialize(HierarchicUnserializer& unserializer);
//     virtual void rescale(int dxScale, int dtScale);
//     Array<T,2> const& getPos1() const { return pos1; }
//     Array<T,2> const& getPos2() const { return pos2; }
//     bool endOfCycle() const { return stage==3; }
// public:
//     virtual bool crossedSurface1() const =0;
//     virtual bool crossedSurface2() const =0;
//     virtual T getSurfaceDistance() const =0;
// private:
//     int stage;
//     Array<T,2> pos1, pos2;
// };
//
// template<typename T, template<typename U> class Descriptor>
// class MappingParticleZslice2D : public MappingParticle2D<T,Descriptor> {
// public:
//     MappingParticleZslice2D();
//     MappingParticleZslice2D(plint tag_, Array<T,2> const& position_, Array<T,2> const& velocity_,
//     T zSlice1_, T zSlice2_); virtual int getId() const; virtual void
//     serialize(HierarchicSerializer& serializer) const; virtual void
//     unserialize(HierarchicUnserializer& unserializer); virtual
//     MappingParticleZslice2D<T,Descriptor>* clone() const; virtual void rescale(int dxScale, int
//     dtScale);
// public:
//     virtual bool crossedSurface1() const;
//     virtual bool crossedSurface2() const;
//     virtual T getSurfaceDistance() const;
// private:
//     T zSlice1, zSlice2;
//     static int id;
// };

namespace meta {

template <typename T, template <typename U> class Descriptor>
ParticleRegistration2D<T, Descriptor> &particleRegistration2D();

template <typename T, template <typename U> class Descriptor, class VisualParticle>
class VisualParticleGenerator2D : public ParticleGenerator2D<T, Descriptor> {
    virtual Particle2D<T, Descriptor> *generate(HierarchicUnserializer &unserializer) const
    {
        // tag, position, scalars, vectors.
        plint tag;
        unserializer.readValue(tag);
        Array<T, 2> position;
        unserializer.readValues<T, 2>(position);
        pluint scalarSize;
        unserializer.readValue(scalarSize);
        std::vector<T> scalars(scalarSize);
        unserializer.readValues(scalars);
        pluint vectorSize;
        unserializer.readValue(vectorSize);
        std::vector<Array<T, 2> > vectors(vectorSize);
        unserializer.readValues(vectors);
        return new VisualParticle(tag, position, scalars, vectors);
    }
};

template <typename T, template <typename U> class Descriptor, class VisualParticle>
int registerVisualParticle2D(std::string name)
{
    return particleRegistration2D<T, Descriptor>().announce(
        name, new VisualParticleGenerator2D<T, Descriptor, VisualParticle>);
}

}  // namespace meta

}  // namespace plb

#endif  // VISUAL_PARTICLE_2D_H
