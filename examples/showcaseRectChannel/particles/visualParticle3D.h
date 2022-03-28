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

#ifndef VISUAL_PARTICLE_3D_H
#define VISUAL_PARTICLE_3D_H

#include <vector>

#include "atomicBlock/blockLattice3D.h"
#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "particles/particle3D.h"
#include "particles/particleIdentifiers3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class VisualParticle3D : public Particle3D<T, Descriptor> {
public:
    VisualParticle3D();
    VisualParticle3D(plint tag_, Array<T, 3> const &position);
    VisualParticle3D(
        plint tag_, Array<T, 3> const &position, std::vector<T> const &scalars_,
        std::vector<Array<T, 3> > const &vectors_);
    virtual void velocityToParticle(TensorField3D<T, 3> &velocityField, T scaling = 1.) { }
    virtual void velocityToParticle(NTensorField3D<T> &velocityField, T scaling = 1.) { }
    virtual void rhoBarJtoParticle(NTensorField3D<T> &rhoBarJfield, bool velIsJ, T scaling = 1.) { }
    virtual void fluidToParticle(BlockLattice3D<T, Descriptor> &fluid, T scaling = 1.) { }
    virtual void fluidToParticle(
        BlockLattice3D<T, Descriptor> &fluid, ScalarField3D<T> &sF, T scaling = 1.)
    { }
    virtual void advance() { }
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual int getId() const;
    virtual void reset(Array<T, 3> const &position);
    virtual VisualParticle3D<T, Descriptor> *clone() const;
    virtual bool getScalar(plint whichScalar, T &scalar) const;
    virtual bool setScalar(plint whichScalar, T scalar);
    virtual bool setScalars(std::vector<T> const &scalars_);
    virtual bool getVector(plint whichVector, Array<T, 3> &vector) const;
    virtual bool setVector(plint whichVector, Array<T, 3> const &vector);
    virtual bool setVectors(std::vector<Array<T, 3> > const &vectors_);

private:
    std::vector<T> scalars;
    std::vector<Array<T, 3> > vectors;
    static int id;
};

template <typename T, template <typename U> class Descriptor>
class MappingParticle3D : public PointParticle3D<T, Descriptor> {
public:
    MappingParticle3D();
    MappingParticle3D(plint tag_, Array<T, 3> const &position_, Array<T, 3> const &velocity_);
    virtual void advance();
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual void rescale(int dxScale, int dtScale);
    Array<T, 3> const &getPos1() const
    {
        return pos1;
    }
    Array<T, 3> const &getPos2() const
    {
        return pos2;
    }
    bool endOfCycle() const
    {
        return stage == 3;
    }

public:
    virtual bool crossedSurface1() const = 0;
    virtual bool crossedSurface2() const = 0;
    virtual T getSurfaceDistance() const = 0;

private:
    int stage;
    Array<T, 3> pos1, pos2;
};

template <typename T, template <typename U> class Descriptor>
class MappingParticleZslice3D : public MappingParticle3D<T, Descriptor> {
public:
    MappingParticleZslice3D();
    MappingParticleZslice3D(
        plint tag_, Array<T, 3> const &position_, Array<T, 3> const &velocity_, T zSlice1_,
        T zSlice2_);
    virtual int getId() const;
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual MappingParticleZslice3D<T, Descriptor> *clone() const;
    virtual void rescale(int dxScale, int dtScale);

public:
    virtual bool crossedSurface1() const;
    virtual bool crossedSurface2() const;
    virtual T getSurfaceDistance() const;

private:
    T zSlice1, zSlice2;
    static int id;
};

/* PlaneMappingParticle3D: This particle type, is used essentially to compute
 * streamlines between an initial particle position, and a position where the
 * particle meets a terminal plane. It also keeps the residence time, which is
 * the difference between the current time and the time the particle was injected.
 * When the particle reaches the terminal plane, its position and residence time
 * freeze. The time scaling variable serves as the ratio between the time step
 * of the particle integration, and the time step of the fluid integration.
 * */
template <typename T, template <typename U> class Descriptor>
class PlaneMappingParticle3D : public PointParticle3D<T, Descriptor> {
public:
    PlaneMappingParticle3D();
    PlaneMappingParticle3D(
        plint tag_, Array<T, 3> const &position_, Array<T, 3> const &velocity_,
        Plane<T> const &terminalPlane_, bool advanceBackwardInTime_, T timeScaling_);
    virtual PlaneMappingParticle3D<T, Descriptor> *clone() const;
    virtual int getId() const;
    virtual void advance();
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual void rescale(int dxScale, int dtScale);
    Array<T, 3> const &getInitialPosition() const
    {
        return initialPosition;
    }
    Array<T, 3> const &getInitialVelocity() const
    {
        return initialVelocity;
    }
    Array<T, 3> &getInitialVelocity()
    {
        return initialVelocity;
    }
    Array<T, 3> const &getTerminalVelocity() const
    {
        return terminalVelocity;
    }
    Plane<T> const &getTerminalPlane() const
    {
        return terminalPlane;
    }
    T getResidenceTime() const
    {
        return residenceTime;
    }
    T getTimeScaling() const
    {
        return timeScaling;
    }
    T &getTimeScaling()
    {
        return timeScaling;
    }
    virtual bool passedTerminalPlane() const
    {
        return reachedTerminalPlane;
    }

private:
    Plane<T> terminalPlane;
    Array<T, 3> initialPosition;
    Array<T, 3> initialVelocity;
    Array<T, 3> terminalVelocity;
    T residenceTime;
    bool advanceBackwardInTime;
    T timeScaling;
    bool reachedTerminalPlane;
    static int id;
};

/*
 * TimeRegisteringParticle3D: This particle type, stores the time it exists
 * (its age) since passing from a user-defined plane. The way this particle is
 * supposed to be used is the following. In the beginning the particle is
 * injected at a position which belongs to the half space A. The time registering
 * has not started yet. As the particle moves, at some point it will cross the
 * plane and will pass to the half space B. Then the time registering begins, and
 * the crossing of the plane is never again checked. The normal of the plane must
 * point from half space A to B. The time scaling variable serves as the ratio
 * between the time step of the particle integration, and the time step of the
 * fluid integration.
 */
template <typename T, template <typename U> class Descriptor>
class TimeRegisteringParticle3D : public PointParticle3D<T, Descriptor> {
public:
    TimeRegisteringParticle3D();
    TimeRegisteringParticle3D(
        plint tag_, Array<T, 3> const &position_, Array<T, 3> const &velocity_,
        Plane<T> const &initialPlane_, T timeScaling_);
    virtual TimeRegisteringParticle3D<T, Descriptor> *clone() const;
    virtual int getId() const;
    virtual void advance();
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    virtual void rescale(int dxScale, int dtScale);
    Plane<T> const &getInitialPlane() const
    {
        return initialPlane;
    }
    T getRegisteredTime() const
    {
        return registeredTime;
    }
    T getTimeScaling() const
    {
        return timeScaling;
    }
    T &getTimeScaling()
    {
        return timeScaling;
    }
    virtual bool passedInitialPlane() const
    {
        return reachedInitialPlane;
    }

private:
    Plane<T> initialPlane;
    T registeredTime;
    T timeScaling;
    bool reachedInitialPlane;
    static int id;
};

namespace meta {

template <typename T, template <typename U> class Descriptor>
ParticleRegistration3D<T, Descriptor> &particleRegistration3D();

template <typename T, template <typename U> class Descriptor, class VisualParticle>
class VisualParticleGenerator3D : public ParticleGenerator3D<T, Descriptor> {
    virtual Particle3D<T, Descriptor> *generate(HierarchicUnserializer &unserializer) const
    {
        // tag, position, scalars, vectors.
        plint tag;
        unserializer.readValue(tag);
        Array<T, 3> position;
        unserializer.readValues<T, 3>(position);
        pluint scalarSize;
        unserializer.readValue(scalarSize);
        std::vector<T> scalars(scalarSize);
        unserializer.readValues(scalars);
        pluint vectorSize;
        unserializer.readValue(vectorSize);
        std::vector<Array<T, 3> > vectors(vectorSize);
        unserializer.readValues(vectors);
        return new VisualParticle(tag, position, scalars, vectors);
    }
};

template <typename T, template <typename U> class Descriptor, class VisualParticle>
int registerVisualParticle3D(std::string name)
{
    return particleRegistration3D<T, Descriptor>().announce(
        name, new VisualParticleGenerator3D<T, Descriptor, VisualParticle>);
}

}  // namespace meta

}  // namespace plb

#endif  // VISUAL_PARTICLE_3D_H
