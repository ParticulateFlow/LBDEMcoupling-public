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

#ifndef PARTICLE_3D_HH
#define PARTICLE_3D_HH

#include <cmath>

#include "core/globalDefs.h"
#include "finiteDifference/interpolations3D.h"
#include "multiGrid/multiGridUtil.h"
#include "particles/particle3D.h"
#include "particles/particleIdentifiers3D.h"

namespace plb {

/* *************** class Particle3D ***************************************** */

template <typename T, template <typename U> class Descriptor>
Particle3D<T, Descriptor>::Particle3D() : tag(0), position(T(), T(), T())
{ }

template <typename T, template <typename U> class Descriptor>
Particle3D<T, Descriptor>::Particle3D(plint tag_, Array<T, 3> const &position_) :
    tag(tag_), position(position_)
{ }

template <typename T, template <typename U> class Descriptor>
void Particle3D<T, Descriptor>::reset(Array<T, 3> const &position_)
{
    position = position_;
}

template <typename T, template <typename U> class Descriptor>
void Particle3D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    serializer.addValue(tag);
    serializer.addValues<T, 3>(position);
}

template <typename T, template <typename U> class Descriptor>
void Particle3D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    unserializer.readValue(tag);
    unserializer.readValues<T, 3>(position);
}

template <typename T, template <typename U> class Descriptor>
plint Particle3D<T, Descriptor>::getTag() const
{
    return tag;
}

template <typename T, template <typename U> class Descriptor>
void Particle3D<T, Descriptor>::setTag(plint tag_)
{
    tag = tag_;
}

template <typename T, template <typename U> class Descriptor>
bool Particle3D<T, Descriptor>::getVector(plint whichVector, Array<T, 3> &vector) const
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle3D<T, Descriptor>::setVector(plint whichVector, Array<T, 3> const &vector)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle3D<T, Descriptor>::getScalar(plint whichScalar, T &scalar) const
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle3D<T, Descriptor>::setScalar(plint whichScalar, T scalar)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle3D<T, Descriptor>::getTensor(
    plint whichVector, Array<T, SymmetricTensorImpl<T, 3>::n> &tensor) const
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle3D<T, Descriptor>::setTensor(
    plint whichVector, Array<T, SymmetricTensorImpl<T, 3>::n> const &tensor)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle3D<T, Descriptor>::setScalars(std::vector<T> const &scalars)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle3D<T, Descriptor>::setVectors(std::vector<Array<T, 3> > const &vectors)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle3D<T, Descriptor>::setTensors(
    std::vector<Array<T, SymmetricTensorImpl<T, 3>::n> > const &tensors)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
void Particle3D<T, Descriptor>::rescale(int dxScale, int dtScale)
{
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    position *= scaleFactor;
}

/* *************** class PointParticle3D ************************************ */

template <typename T, template <typename U> class Descriptor>
int PointParticle3D<T, Descriptor>::id =
    meta::registerPointParticle3D<T, Descriptor, PointParticle3D<T, Descriptor> >("Point");

template <typename T, template <typename U> class Descriptor>
PointParticle3D<T, Descriptor>::PointParticle3D() :
    Particle3D<T, Descriptor>(), velocity(T(), T(), T())
{ }

template <typename T, template <typename U> class Descriptor>
PointParticle3D<T, Descriptor>::PointParticle3D(
    plint tag_, Array<T, 3> const &position_, Array<T, 3> const &velocity_) :
    Particle3D<T, Descriptor>(tag_, position_), velocity(velocity_)
{ }

template <typename T, template <typename U> class Descriptor>
int PointParticle3D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void PointParticle3D<T, Descriptor>::reset(Array<T, 3> const &position_)
{
    Particle3D<T, Descriptor>::reset(position_);
    velocity.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
void PointParticle3D<T, Descriptor>::velocityToParticle(
    TensorField3D<T, 3> &velocityField, T scaling)
{
    velocity = predictorCorrectorTensorField<T, 3>(velocityField, this->getPosition(), scaling);
}

template <typename T, template <typename U> class Descriptor>
void PointParticle3D<T, Descriptor>::velocityToParticle(NTensorField3D<T> &velocityField, T scaling)
{
    velocity = predictorCorrectorNTensorField<T>(velocityField, this->getPosition(), scaling);
}

template <typename T, template <typename U> class Descriptor>
void PointParticle3D<T, Descriptor>::rhoBarJtoParticle(
    NTensorField3D<T> &rhoBarJfield, bool velIsJ, T scaling)
{
    T rhoBar;
    Array<T, 3> j;
    predictorCorrectorRhoBarJ(rhoBarJfield, this->getPosition(), velIsJ, j, rhoBar);
    if (velIsJ) {
        velocity = j * scaling;
    } else {
        velocity = j * scaling * Descriptor<T>::invRho(rhoBar);
    }
}

template <typename T, template <typename U> class Descriptor>
void PointParticle3D<T, Descriptor>::fluidToParticle(
    BlockLattice3D<T, Descriptor> &fluid, T scaling)
{
    static const T maxVel = 0.25 - 1.e-6;
    static const T maxVelSqr = maxVel * maxVel;
#ifdef PLB_DEBUG
    Box3D bbox(fluid.getBoundingBox());
#endif
    Dot3D loc(fluid.getLocation());
    Array<T, 3> position1(this->getPosition() - Array<T, 3>(loc.x, loc.y, loc.z));
    PLB_ASSERT(position1[0] >= bbox.x0 + 0.5);
    PLB_ASSERT(position1[0] <= bbox.x1 - 0.5);
    PLB_ASSERT(position1[1] >= bbox.y0 + 0.5);
    PLB_ASSERT(position1[1] <= bbox.y1 - 0.5);
    PLB_ASSERT(position1[2] >= bbox.z0 + 0.5);
    PLB_ASSERT(position1[2] <= bbox.z1 - 0.5);

    Dot3D intPos((plint)position1[0], (plint)position1[1], (plint)position1[2]);
    T u = position1[0] - (T)intPos.x;
    T v = position1[1] - (T)intPos.y;
    T w = position1[2] - (T)intPos.z;

    Array<T, 3> tmpVel;
    Array<T, 3> velocity1;
    velocity1.resetToZero();
    fluid.get(intPos.x, intPos.y, intPos.z).computeVelocity(tmpVel);
    velocity1 += ((T)1. - u) * ((T)1. - v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x, intPos.y, intPos.z + 1).computeVelocity(tmpVel);
    velocity1 += ((T)1. - u) * ((T)1. - v) * (w)*tmpVel;
    fluid.get(intPos.x, intPos.y + 1, intPos.z).computeVelocity(tmpVel);
    velocity1 += ((T)1. - u) * (v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x, intPos.y + 1, intPos.z + 1).computeVelocity(tmpVel);
    velocity1 += ((T)1. - u) * (v) * (w)*tmpVel;
    fluid.get(intPos.x + 1, intPos.y, intPos.z).computeVelocity(tmpVel);
    velocity1 += (u) * ((T)1. - v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x + 1, intPos.y, intPos.z + 1).computeVelocity(tmpVel);
    velocity1 += (u) * ((T)1. - v) * (w)*tmpVel;
    fluid.get(intPos.x + 1, intPos.y + 1, intPos.z).computeVelocity(tmpVel);
    velocity1 += (u) * (v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x + 1, intPos.y + 1, intPos.z + 1).computeVelocity(tmpVel);
    velocity1 += (u) * (v) * (w)*tmpVel;

    velocity1 *= scaling;

    if (normSqr(velocity1) > maxVelSqr) {
        velocity1 /= norm(velocity1);
        velocity1 *= maxVel;
    }

    Array<T, 3> position2(position1 + velocity1);

    intPos = Dot3D((plint)position2[0], (plint)position2[1], (plint)position2[2]);
    u = position2[0] - (T)intPos.x;
    v = position2[1] - (T)intPos.y;
    w = position2[2] - (T)intPos.z;

    Array<T, 3> velocity2;
    velocity2.resetToZero();
    fluid.get(intPos.x, intPos.y, intPos.z).computeVelocity(tmpVel);
    velocity2 += ((T)1. - u) * ((T)1. - v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x, intPos.y, intPos.z + 1).computeVelocity(tmpVel);
    velocity2 += ((T)1. - u) * ((T)1. - v) * (w)*tmpVel;
    fluid.get(intPos.x, intPos.y + 1, intPos.z).computeVelocity(tmpVel);
    velocity2 += ((T)1. - u) * (v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x, intPos.y + 1, intPos.z + 1).computeVelocity(tmpVel);
    velocity2 += ((T)1. - u) * (v) * (w)*tmpVel;
    fluid.get(intPos.x + 1, intPos.y, intPos.z).computeVelocity(tmpVel);
    velocity2 += (u) * ((T)1. - v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x + 1, intPos.y, intPos.z + 1).computeVelocity(tmpVel);
    velocity2 += (u) * ((T)1. - v) * (w)*tmpVel;
    fluid.get(intPos.x + 1, intPos.y + 1, intPos.z).computeVelocity(tmpVel);
    velocity2 += (u) * (v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x + 1, intPos.y + 1, intPos.z + 1).computeVelocity(tmpVel);
    velocity2 += (u) * (v) * (w)*tmpVel;

    velocity2 *= scaling;

    if (normSqr(velocity2) > maxVelSqr) {
        velocity2 /= norm(velocity2);
        velocity2 *= maxVel;
    }

    velocity = (velocity1 + velocity2) / (T)2;
}

template <typename T, template <typename U> class Descriptor>
void PointParticle3D<T, Descriptor>::fluidToParticle(
    BlockLattice3D<T, Descriptor> &fluid, ScalarField3D<T> &sF, T scaling)
{
    static const T maxVel = 0.25 - 1.e-6;
    static const T maxVelSqr = maxVel * maxVel;
#ifdef PLB_DEBUG
    Box3D bbox(fluid.getBoundingBox());
#endif
    Dot3D loc(fluid.getLocation());
    Array<T, 3> position1(this->getPosition() - Array<T, 3>(loc.x, loc.y, loc.z));
    PLB_ASSERT(position1[0] >= bbox.x0 + 0.5);
    PLB_ASSERT(position1[0] <= bbox.x1 - 0.5);
    PLB_ASSERT(position1[1] >= bbox.y0 + 0.5);
    PLB_ASSERT(position1[1] <= bbox.y1 - 0.5);
    PLB_ASSERT(position1[2] >= bbox.z0 + 0.5);
    PLB_ASSERT(position1[2] <= bbox.z1 - 0.5);

    Dot3D intPos((plint)position1[0], (plint)position1[1], (plint)position1[2]);
    T u = position1[0] - (T)intPos.x;
    T v = position1[1] - (T)intPos.y;
    T w = position1[2] - (T)intPos.z;

    Array<T, 3> tmpVel;
    Array<T, 3> velocity1;
    velocity1.resetToZero();
    fluid.get(intPos.x, intPos.y, intPos.z).computeVelocity(tmpVel);
    velocity1 += ((T)1. - u) * ((T)1. - v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x, intPos.y, intPos.z + 1).computeVelocity(tmpVel);
    velocity1 += ((T)1. - u) * ((T)1. - v) * (w)*tmpVel;
    fluid.get(intPos.x, intPos.y + 1, intPos.z).computeVelocity(tmpVel);
    velocity1 += ((T)1. - u) * (v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x, intPos.y + 1, intPos.z + 1).computeVelocity(tmpVel);
    velocity1 += ((T)1. - u) * (v) * (w)*tmpVel;
    fluid.get(intPos.x + 1, intPos.y, intPos.z).computeVelocity(tmpVel);
    velocity1 += (u) * ((T)1. - v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x + 1, intPos.y, intPos.z + 1).computeVelocity(tmpVel);
    velocity1 += (u) * ((T)1. - v) * (w)*tmpVel;
    fluid.get(intPos.x + 1, intPos.y + 1, intPos.z).computeVelocity(tmpVel);
    velocity1 += (u) * (v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x + 1, intPos.y + 1, intPos.z + 1).computeVelocity(tmpVel);
    velocity1 += (u) * (v) * (w)*tmpVel;

    velocity1 *= scaling;

    if (normSqr(velocity1) > maxVelSqr) {
        velocity1 /= norm(velocity1);
        velocity1 *= maxVel;
    }

    Array<T, 3> position2(position1 + velocity1);

    intPos = Dot3D((plint)position2[0], (plint)position2[1], (plint)position2[2]);
    u = position2[0] - (T)intPos.x;
    v = position2[1] - (T)intPos.y;
    w = position2[2] - (T)intPos.z;

    Array<T, 3> velocity2;
    velocity2.resetToZero();
    fluid.get(intPos.x, intPos.y, intPos.z).computeVelocity(tmpVel);
    velocity2 += ((T)1. - u) * ((T)1. - v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x, intPos.y, intPos.z + 1).computeVelocity(tmpVel);
    velocity2 += ((T)1. - u) * ((T)1. - v) * (w)*tmpVel;
    fluid.get(intPos.x, intPos.y + 1, intPos.z).computeVelocity(tmpVel);
    velocity2 += ((T)1. - u) * (v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x, intPos.y + 1, intPos.z + 1).computeVelocity(tmpVel);
    velocity2 += ((T)1. - u) * (v) * (w)*tmpVel;
    fluid.get(intPos.x + 1, intPos.y, intPos.z).computeVelocity(tmpVel);
    velocity2 += (u) * ((T)1. - v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x + 1, intPos.y, intPos.z + 1).computeVelocity(tmpVel);
    velocity2 += (u) * ((T)1. - v) * (w)*tmpVel;
    fluid.get(intPos.x + 1, intPos.y + 1, intPos.z).computeVelocity(tmpVel);
    velocity2 += (u) * (v) * ((T)1. - w) * tmpVel;
    fluid.get(intPos.x + 1, intPos.y + 1, intPos.z + 1).computeVelocity(tmpVel);
    velocity2 += (u) * (v) * (w)*tmpVel;

    velocity2 *= scaling;

    if (normSqr(velocity2) > maxVelSqr) {
        velocity2 /= norm(velocity2);
        velocity2 *= maxVel;
    }
    if (contained(intPos, sF.getBoundingBox()))
        velocity = (velocity1 + velocity2) * (1. - sF.get(intPos.x, intPos.y, intPos.z)) / (T)2;

    else
        velocity = (velocity1 + velocity2) / (T)2;
}

template <typename T, template <typename U> class Descriptor>
void PointParticle3D<T, Descriptor>::advance()
{
    PLB_ASSERT(norm(velocity) < 1.);
    this->getPosition() += velocity;
}

template <typename T, template <typename U> class Descriptor>
void PointParticle3D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    Particle3D<T, Descriptor>::serialize(serializer);
    serializer.addValues<T, 3>(velocity);
}

template <typename T, template <typename U> class Descriptor>
void PointParticle3D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    Particle3D<T, Descriptor>::unserialize(unserializer);
    unserializer.readValues<T, 3>(velocity);
}

template <typename T, template <typename U> class Descriptor>
PointParticle3D<T, Descriptor> *PointParticle3D<T, Descriptor>::clone() const
{
    return new PointParticle3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
bool PointParticle3D<T, Descriptor>::getVector(plint whichVector, Array<T, 3> &vector) const
{
    if (whichVector == 0) {
        vector = velocity;
        return true;
    }
    return Particle3D<T, Descriptor>::getVector(whichVector, vector);
}

template <typename T, template <typename U> class Descriptor>
bool PointParticle3D<T, Descriptor>::setVectors(std::vector<Array<T, 3> > const &vectors)
{
    if (vectors.size() == 1) {
        velocity = vectors[0];
        return true;
    }
    return false;
}

template <typename T, template <typename U> class Descriptor>
void PointParticle3D<T, Descriptor>::rescale(int dxScale, int dtScale)
{
    Particle3D<T, Descriptor>::rescale(dxScale, dtScale);
    T dimDx = 1.;
    T dimDt = -1.;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    velocity *= scaleFactor;
}

/* *************** class NormedVelocityParticle3D ************************************ */

template <typename T, template <typename U> class Descriptor>
int NormedVelocityParticle3D<T, Descriptor>::id =
    meta::registerGenericParticle3D<T, Descriptor, NormedVelocityParticle3D<T, Descriptor> >(
        "NormedVelocity");

template <typename T, template <typename U> class Descriptor>
NormedVelocityParticle3D<T, Descriptor>::NormedVelocityParticle3D() :
    PointParticle3D<T, Descriptor>(), fluidUmax(0.), particleUmax(0.), exponent(0.)
{ }

template <typename T, template <typename U> class Descriptor>
NormedVelocityParticle3D<T, Descriptor>::NormedVelocityParticle3D(
    plint tag_, Array<T, 3> const &position_, Array<T, 3> const &velocity_, T fluidUmax_,
    T particleUmax_, T exponent_) :
    PointParticle3D<T, Descriptor>(tag_, position_, velocity_),
    fluidUmax(fluidUmax_),
    particleUmax(particleUmax_),
    exponent(exponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void NormedVelocityParticle3D<T, Descriptor>::advance()
{
    T fluidNormu = norm(this->getVelocity());
    PLB_ASSERT(fluidNormu < 1.);
    T particleNormU;
    if (fluidNormu >= fluidUmax) {
        particleNormU = particleUmax;
    } else {
        particleNormU = std::pow(fluidNormu / fluidUmax, (T)exponent) * particleUmax;
    }
    PLB_ASSERT(particleNormU < 1.);
    this->getPosition() += this->getVelocity() / fluidNormu * particleNormU;
}

template <typename T, template <typename U> class Descriptor>
void NormedVelocityParticle3D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    PointParticle3D<T, Descriptor>::serialize(serializer);
    serializer.addValue(fluidUmax);
    serializer.addValue(particleUmax);
    serializer.addValue(exponent);
}

template <typename T, template <typename U> class Descriptor>
void NormedVelocityParticle3D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    PointParticle3D<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(fluidUmax);
    unserializer.readValue(particleUmax);
    unserializer.readValue(exponent);
}

template <typename T, template <typename U> class Descriptor>
NormedVelocityParticle3D<T, Descriptor> *NormedVelocityParticle3D<T, Descriptor>::clone() const
{
    return new NormedVelocityParticle3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void NormedVelocityParticle3D<T, Descriptor>::rescale(int dxScale, int dtScale)
{
    PointParticle3D<T, Descriptor>::rescale(dxScale, dtScale);
    T dimDx = 1.;
    T dimDt = -1.;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    fluidUmax *= scaleFactor;
    particleUmax *= scaleFactor;
}

template <typename T, template <typename U> class Descriptor>
int NormedVelocityParticle3D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void serialize(Particle3D<T, Descriptor> const &particle, std::vector<char> &data)
{
    HierarchicSerializer serializer(data, particle.getId());
    particle.serialize(serializer);
}

template <typename T, template <typename U> class Descriptor>
void generateAndUnserializeParticles(
    std::vector<char> const &data, std::vector<Particle3D<T, Descriptor> *> &particles)
{
    pluint serializerPos = 0;
    while (serializerPos < data.size()) {
        HierarchicUnserializer unserializer(data, serializerPos);
        particles.push_back(meta::particleRegistration3D<T, Descriptor>().generate(unserializer));
        serializerPos = unserializer.getCurrentPos();
    }
}

/* *************** class RestParticle3D ************************************ */

template <typename T, template <typename U> class Descriptor>
int RestParticle3D<T, Descriptor>::id =
    meta::registerGenericParticle3D<T, Descriptor, RestParticle3D<T, Descriptor> >("Rest");

template <typename T, template <typename U> class Descriptor>
RestParticle3D<T, Descriptor>::RestParticle3D()
{ }

template <typename T, template <typename U> class Descriptor>
RestParticle3D<T, Descriptor>::RestParticle3D(plint tag_, Array<T, 3> const &position) :
    Particle3D<T, Descriptor>(tag_, position)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle3D<T, Descriptor>::velocityToParticle(
    TensorField3D<T, 3> &velocityField, T scaling)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle3D<T, Descriptor>::velocityToParticle(NTensorField3D<T> &velocityField, T scaling)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle3D<T, Descriptor>::rhoBarJtoParticle(
    NTensorField3D<T> &rhoBarJfield, bool velIsJ, T scaling)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle3D<T, Descriptor>::fluidToParticle(BlockLattice3D<T, Descriptor> &fluid, T scaling)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle3D<T, Descriptor>::fluidToParticle(
    BlockLattice3D<T, Descriptor> &fluid, ScalarField3D<T> &sF, T scaling)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle3D<T, Descriptor>::advance()
{ }

template <typename T, template <typename U> class Descriptor>
int RestParticle3D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
RestParticle3D<T, Descriptor> *RestParticle3D<T, Descriptor>::clone() const
{
    return new RestParticle3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
bool RestParticle3D<T, Descriptor>::getVector(plint whichVector, Array<T, 3> &vector) const
{
    return Particle3D<T, Descriptor>::getVector(whichVector, vector);
}

/* *************** class VerletParticle3D ************************************ */

template <typename T, template <typename U> class Descriptor>
int VerletParticle3D<T, Descriptor>::id =
    meta::registerGenericParticle3D<T, Descriptor, VerletParticle3D<T, Descriptor> >("Verlet");

template <typename T, template <typename U> class Descriptor>
VerletParticle3D<T, Descriptor>::VerletParticle3D() :
    v(T(), T(), T()),
    vHalfTime(T(), T(), T()),
    a(T(), T(), T()),
    fluidCompliance((T)1.),
    rho((T)1.),
    invRho((T)1.)
{ }

template <typename T, template <typename U> class Descriptor>
VerletParticle3D<T, Descriptor>::VerletParticle3D(plint tag_, Array<T, 3> const &position) :
    Particle3D<T, Descriptor>(tag_, position),
    v(T(), T(), T()),
    vHalfTime(T(), T(), T()),
    a(T(), T(), T()),
    fluidCompliance((T)1.),
    rho((T)1.),
    invRho((T)1.)
{ }

template <typename T, template <typename U> class Descriptor>
void VerletParticle3D<T, Descriptor>::velocityToParticle(
    TensorField3D<T, 3> &velocityField, T scaling)
{
    Array<T, 3> position(this->getPosition());
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(velocityField, position, pos, weights);
    Array<T, 3> fluidVelocity;
    fluidVelocity.resetToZero();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        fluidVelocity +=
            weights[iCell] * velocityField.get(pos[iCell].x, pos[iCell].y, pos[iCell].z) * scaling;
    }

    Array<T, 3> force((fluidVelocity - this->get_v()) * fluidCompliance);
    this->set_a(force / this->get_rho());
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle3D<T, Descriptor>::velocityToParticle(
    NTensorField3D<T> &velocityField, T scaling)
{
    Array<T, 3> position(this->getPosition());
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(velocityField, position, pos, weights);
    Array<T, 3> fluidVelocity;
    fluidVelocity.resetToZero();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        T *data = velocityField.get(pos[iCell].x, pos[iCell].y, pos[iCell].z);
        fluidVelocity[0] += weights[iCell] * scaling * data[0];
        fluidVelocity[1] += weights[iCell] * scaling * data[1];
        fluidVelocity[2] += weights[iCell] * scaling * data[2];
    }

    Array<T, 3> force((fluidVelocity - this->get_v()) * fluidCompliance);
    this->set_a(force / this->get_rho());
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle3D<T, Descriptor>::rhoBarJtoParticle(
    NTensorField3D<T> &rhoBarJfield, bool velIsJ, T scaling)
{
    Array<T, 3> position(this->getPosition());
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(rhoBarJfield, position, pos, weights);
    Array<T, 3> j;
    j.resetToZero();
    T rhoBar = T();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        T *data = rhoBarJfield.get(pos[iCell].x, pos[iCell].y, pos[iCell].z);
        j.add_from_cArray(data + 1, weights[iCell]);
        rhoBar += weights[iCell] * (*data);
    }
    Array<T, 3> fluidVelocity;
    if (velIsJ) {
        fluidVelocity = j * scaling;
    } else {
        fluidVelocity = j * scaling * Descriptor<T>::invRho(rhoBar);
    }

    Array<T, 3> force((fluidVelocity - this->get_v()) * fluidCompliance);
    this->set_a(force / this->get_rho());
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle3D<T, Descriptor>::fluidToParticle(
    BlockLattice3D<T, Descriptor> &fluid, T scaling)
{
    Array<T, 3> position(this->getPosition());
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(fluid, position, pos, weights);
    Array<T, 3> tmpVel;
    Array<T, 3> fluidVelocity;
    fluidVelocity.resetToZero();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        if (contained(pos[iCell].x, pos[iCell].y, pos[iCell].z, fluid.getBoundingBox())) {
            fluid.get(pos[iCell].x, pos[iCell].y, pos[iCell].z).computeVelocity(tmpVel);
        } else {
            tmpVel.resetToZero();
        }
        fluidVelocity += weights[iCell] * tmpVel * scaling;
    }

    Array<T, 3> force((fluidVelocity - this->get_v()) * fluidCompliance);
    this->set_a(force / this->get_rho());
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle3D<T, Descriptor>::fluidToParticle(
    BlockLattice3D<T, Descriptor> &fluid, ScalarField3D<T> &sF, T scaling)
{
    Array<T, 3> position(this->getPosition());
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(fluid, position, pos, weights);
    Array<T, 3> tmpVel;
    Array<T, 3> fluidVelocity;
    fluidVelocity.resetToZero();
    Dot3D intPos((plint)position[0], (plint)position[1], (plint)position[2]);
    for (plint iCell = 0; iCell < 8; ++iCell) {
        if (contained(pos[iCell].x, pos[iCell].y, pos[iCell].z, fluid.getBoundingBox())) {
            fluid.get(pos[iCell].x, pos[iCell].y, pos[iCell].z).computeVelocity(tmpVel);
        } else {
            tmpVel.resetToZero();
        }
        if (contained(intPos, sF.getBoundingBox()))
            fluidVelocity +=
                weights[iCell] * tmpVel * scaling * (1. - sF.get(intPos.x, intPos.y, intPos.z));
        else
            fluidVelocity += weights[iCell] * tmpVel * scaling;
    }

    Array<T, 3> force((fluidVelocity - this->get_v()) * fluidCompliance);
    this->set_a(force / this->get_rho());
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle3D<T, Descriptor>::advance()
{
    vHalfTime = v + (T)0.5 * a;
    this->getPosition() += vHalfTime;
}

template <typename T, template <typename U> class Descriptor>
int VerletParticle3D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle3D<T, Descriptor>::reset(Array<T, 3> const &position_)
{
    Particle3D<T, Descriptor>::reset(position_);
    v.resetToZero();
    vHalfTime.resetToZero();
    a.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle3D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    Particle3D<T, Descriptor>::serialize(serializer);
    serializer.addValues<T, 3>(v);
    serializer.addValues<T, 3>(vHalfTime);
    serializer.addValues<T, 3>(a);
    serializer.addValue<T>(fluidCompliance);
    serializer.addValue<T>(rho);
    serializer.addValue<T>(invRho);
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle3D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    Particle3D<T, Descriptor>::unserialize(unserializer);
    unserializer.readValues<T, 3>(v);
    unserializer.readValues<T, 3>(vHalfTime);
    unserializer.readValues<T, 3>(a);
    unserializer.readValue<T>(fluidCompliance);
    unserializer.readValue<T>(rho);
    unserializer.readValue<T>(invRho);
}

template <typename T, template <typename U> class Descriptor>
VerletParticle3D<T, Descriptor> *VerletParticle3D<T, Descriptor>::clone() const
{
    return new VerletParticle3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
bool VerletParticle3D<T, Descriptor>::getVector(plint whichVector, Array<T, 3> &vector) const
{
    if (whichVector == 0) {
        vector = get_v();
        return true;
    } else if (whichVector == 1) {
        vector = get_a();
        return true;
    }
    return Particle3D<T, Descriptor>::getVector(whichVector, vector);
}

template <typename T, template <typename U> class Descriptor>
bool VerletParticle3D<T, Descriptor>::getScalar(plint whichScalar, T &scalar) const
{
    if (whichScalar == 0) {
        scalar = rho;
        return true;
    } else if (whichScalar == 1) {
        scalar = fluidCompliance;
        return true;
    }
    return Particle3D<T, Descriptor>::getScalar(whichScalar, scalar);
}

}  // namespace plb

#endif  // PARTICLE_3D_H
