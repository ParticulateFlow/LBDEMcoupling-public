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

#ifndef PARTICLE_2D_HH
#define PARTICLE_2D_HH

#include <cmath>

#include "core/globalDefs.h"
#include "finiteDifference/interpolations2D.h"
#include "particles/particle2D.h"
#include "particles/particleIdentifiers2D.h"

namespace plb {

/* *************** class Particle2D ***************************************** */

template <typename T, template <typename U> class Descriptor>
Particle2D<T, Descriptor>::Particle2D() : tag(0), position(T(), T())
{ }

template <typename T, template <typename U> class Descriptor>
Particle2D<T, Descriptor>::Particle2D(plint tag_, Array<T, 2> const &position_) :
    tag(tag_), position(position_)
{ }

template <typename T, template <typename U> class Descriptor>
void Particle2D<T, Descriptor>::reset(Array<T, 2> const &position_)
{
    position = position_;
}

template <typename T, template <typename U> class Descriptor>
void Particle2D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    serializer.addValue(tag);
    serializer.addValues<T, 2>(position);
}

template <typename T, template <typename U> class Descriptor>
void Particle2D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    unserializer.readValue(tag);
    unserializer.readValues<T, 2>(position);
}

template <typename T, template <typename U> class Descriptor>
plint Particle2D<T, Descriptor>::getTag() const
{
    return tag;
}

template <typename T, template <typename U> class Descriptor>
void Particle2D<T, Descriptor>::setTag(plint tag_)
{
    tag = tag_;
}

template <typename T, template <typename U> class Descriptor>
bool Particle2D<T, Descriptor>::getVector(plint whichVector, Array<T, 2> &vector) const
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle2D<T, Descriptor>::setVector(plint whichVector, Array<T, 2> const &vector)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle2D<T, Descriptor>::getScalar(plint whichScalar, T &scalar) const
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle2D<T, Descriptor>::setScalar(plint whichScalar, T scalar)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle2D<T, Descriptor>::getTensor(
    plint whichVector, Array<T, SymmetricTensorImpl<T, 2>::n> &tensor) const
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle2D<T, Descriptor>::setTensor(
    plint whichVector, Array<T, SymmetricTensorImpl<T, 2>::n> const &tensor)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle2D<T, Descriptor>::setScalars(std::vector<T> const &scalars)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle2D<T, Descriptor>::setVectors(std::vector<Array<T, 2> > const &vectors)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
bool Particle2D<T, Descriptor>::setTensors(
    std::vector<Array<T, SymmetricTensorImpl<T, 2>::n> > const &tensors)
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
void Particle2D<T, Descriptor>::rescale(int dxScale, int dtScale)
{
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    position *= scaleFactor;
}

/* *************** class PointParticle2D ************************************ */

template <typename T, template <typename U> class Descriptor>
int PointParticle2D<T, Descriptor>::id =
    meta::registerPointParticle2D<T, Descriptor, PointParticle2D<T, Descriptor> >("Point");

template <typename T, template <typename U> class Descriptor>
PointParticle2D<T, Descriptor>::PointParticle2D() : Particle2D<T, Descriptor>(), velocity(T(), T())
{ }

template <typename T, template <typename U> class Descriptor>
PointParticle2D<T, Descriptor>::PointParticle2D(
    plint tag_, Array<T, 2> const &position_, Array<T, 2> const &velocity_) :
    Particle2D<T, Descriptor>(tag_, position_), velocity(velocity_)
{ }

template <typename T, template <typename U> class Descriptor>
int PointParticle2D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void PointParticle2D<T, Descriptor>::reset(Array<T, 2> const &position_)
{
    Particle2D<T, Descriptor>::reset(position_);
    velocity.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
void PointParticle2D<T, Descriptor>::velocityToParticle(
    TensorField2D<T, 2> &velocityField, T scaling)
{
    velocity = predictorCorrectorTensorField<T, 2>(velocityField, this->getPosition(), scaling);
}

template <typename T, template <typename U> class Descriptor>
void PointParticle2D<T, Descriptor>::velocityToParticle(NTensorField2D<T> &velocityField, T scaling)
{
    velocity = predictorCorrectorNTensorField<T>(velocityField, this->getPosition(), scaling);
}

template <typename T, template <typename U> class Descriptor>
void PointParticle2D<T, Descriptor>::rhoBarJtoParticle(
    NTensorField2D<T> &rhoBarJfield, bool velIsJ, T scaling)
{
    T rhoBar;
    Array<T, 2> j;
    predictorCorrectorRhoBarJ(rhoBarJfield, this->getPosition(), velIsJ, j, rhoBar);
    if (velIsJ) {
        velocity = j * scaling;
    } else {
        velocity = j * scaling * Descriptor<T>::invRho(rhoBar);
    }
}

template <typename T, template <typename U> class Descriptor>
void PointParticle2D<T, Descriptor>::fluidToParticle(
    BlockLattice2D<T, Descriptor> &fluid, T scaling)
{
    static const T maxVel = 0.25 - 1.e-6;
    static const T maxVelSqr = maxVel * maxVel;
#ifdef PLB_DEBUG
    Box2D bbox(fluid.getBoundingBox());
#endif
    Dot2D loc(fluid.getLocation());
    Array<T, 2> position1(this->getPosition() - Array<T, 2>(loc.x, loc.y));
    PLB_ASSERT(position1[0] >= bbox.x0 + 0.5);
    PLB_ASSERT(position1[0] <= bbox.x1 - 0.5);
    PLB_ASSERT(position1[1] >= bbox.y0 + 0.5);
    PLB_ASSERT(position1[1] <= bbox.y1 - 0.5);

    Dot2D intPos((plint)position1[0], (plint)position1[1]);
    T u = position1[0] - (T)intPos.x;
    T v = position1[1] - (T)intPos.y;

    Array<T, 2> tmpVel;
    Array<T, 2> velocity1;
    velocity1.resetToZero();
    fluid.get(intPos.x, intPos.y).computeVelocity(tmpVel);
    velocity1 += ((T)1. - u) * ((T)1. - v) * tmpVel;
    fluid.get(intPos.x, intPos.y + 1).computeVelocity(tmpVel);
    velocity1 += ((T)1. - u) * (v)*tmpVel;
    fluid.get(intPos.x + 1, intPos.y).computeVelocity(tmpVel);
    velocity1 += (u) * ((T)1. - v) * tmpVel;
    fluid.get(intPos.x + 1, intPos.y + 1).computeVelocity(tmpVel);
    velocity1 += (u) * (v)*tmpVel;

    velocity1 *= scaling;

    if (normSqr(velocity1) > maxVelSqr) {
        velocity1 /= norm(velocity1);
        velocity1 *= maxVel;
    }

    Array<T, 2> position2(position1 + velocity1);

    intPos = Dot2D((plint)position2[0], (plint)position2[1]);
    u = position2[0] - (T)intPos.x;
    v = position2[1] - (T)intPos.y;

    Array<T, 2> velocity2;
    velocity2.resetToZero();
    fluid.get(intPos.x, intPos.y).computeVelocity(tmpVel);
    velocity2 += ((T)1. - u) * ((T)1. - v) * tmpVel;
    fluid.get(intPos.x, intPos.y + 1).computeVelocity(tmpVel);
    velocity2 += ((T)1. - u) * (v)*tmpVel;
    fluid.get(intPos.x + 1, intPos.y).computeVelocity(tmpVel);
    velocity2 += (u) * ((T)1. - v) * tmpVel;
    fluid.get(intPos.x + 1, intPos.y + 1).computeVelocity(tmpVel);
    velocity2 += (u) * (v)*tmpVel;

    velocity2 *= scaling;

    if (normSqr(velocity2) > maxVelSqr) {
        velocity2 /= norm(velocity2);
        velocity2 *= maxVel;
    }

    velocity = (velocity1 + velocity2) / (T)2;
}

template <typename T, template <typename U> class Descriptor>
void PointParticle2D<T, Descriptor>::advance()
{
    PLB_ASSERT(norm(velocity) < 1.);
    this->getPosition() += velocity;
}

template <typename T, template <typename U> class Descriptor>
void PointParticle2D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    Particle2D<T, Descriptor>::serialize(serializer);
    serializer.addValues<T, 2>(velocity);
}

template <typename T, template <typename U> class Descriptor>
void PointParticle2D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    Particle2D<T, Descriptor>::unserialize(unserializer);
    unserializer.readValues<T, 2>(velocity);
}

template <typename T, template <typename U> class Descriptor>
PointParticle2D<T, Descriptor> *PointParticle2D<T, Descriptor>::clone() const
{
    return new PointParticle2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
bool PointParticle2D<T, Descriptor>::getVector(plint whichVector, Array<T, 2> &vector) const
{
    if (whichVector == 0) {
        vector = velocity;
        return true;
    }
    return Particle2D<T, Descriptor>::getVector(whichVector, vector);
}

template <typename T, template <typename U> class Descriptor>
bool PointParticle2D<T, Descriptor>::setVectors(std::vector<Array<T, 2> > const &vectors)
{
    if (vectors.size() == 1) {
        velocity = vectors[0];
        return true;
    }
    return false;
}

template <typename T, template <typename U> class Descriptor>
void PointParticle2D<T, Descriptor>::rescale(int dxScale, int dtScale)
{
    Particle2D<T, Descriptor>::rescale(dxScale, dtScale);
    T dimDx = 1.;
    T dimDt = -1.;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    velocity *= scaleFactor;
}

/* *************** class NormedVelocityParticle2D ************************************ */

template <typename T, template <typename U> class Descriptor>
int NormedVelocityParticle2D<T, Descriptor>::id =
    meta::registerGenericParticle2D<T, Descriptor, NormedVelocityParticle2D<T, Descriptor> >(
        "NormedVelocity");

template <typename T, template <typename U> class Descriptor>
NormedVelocityParticle2D<T, Descriptor>::NormedVelocityParticle2D() :
    PointParticle2D<T, Descriptor>(), fluidUmax(0.), particleUmax(0.), exponent(0.)
{ }

template <typename T, template <typename U> class Descriptor>
NormedVelocityParticle2D<T, Descriptor>::NormedVelocityParticle2D(
    plint tag_, Array<T, 2> const &position_, Array<T, 2> const &velocity_, T fluidUmax_,
    T particleUmax_, T exponent_) :
    PointParticle2D<T, Descriptor>(tag_, position_, velocity_),
    fluidUmax(fluidUmax_),
    particleUmax(particleUmax_),
    exponent(exponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void NormedVelocityParticle2D<T, Descriptor>::advance()
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
void NormedVelocityParticle2D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    PointParticle2D<T, Descriptor>::serialize(serializer);
    serializer.addValue(fluidUmax);
    serializer.addValue(particleUmax);
    serializer.addValue(exponent);
}

template <typename T, template <typename U> class Descriptor>
void NormedVelocityParticle2D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    PointParticle2D<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(fluidUmax);
    unserializer.readValue(particleUmax);
    unserializer.readValue(exponent);
}

template <typename T, template <typename U> class Descriptor>
NormedVelocityParticle2D<T, Descriptor> *NormedVelocityParticle2D<T, Descriptor>::clone() const
{
    return new NormedVelocityParticle2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void NormedVelocityParticle2D<T, Descriptor>::rescale(int dxScale, int dtScale)
{
    PointParticle2D<T, Descriptor>::rescale(dxScale, dtScale);
    T dimDx = 1.;
    T dimDt = -1.;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    fluidUmax *= scaleFactor;
    particleUmax *= scaleFactor;
}

template <typename T, template <typename U> class Descriptor>
int NormedVelocityParticle2D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void serialize(Particle2D<T, Descriptor> const &particle, std::vector<char> &data)
{
    HierarchicSerializer serializer(data, particle.getId());
    particle.serialize(serializer);
}

template <typename T, template <typename U> class Descriptor>
void generateAndUnserializeParticles(
    std::vector<char> const &data, std::vector<Particle2D<T, Descriptor> *> &particles)
{
    pluint serializerPos = 0;
    while (serializerPos < data.size()) {
        HierarchicUnserializer unserializer(data, serializerPos);
        particles.push_back(meta::particleRegistration2D<T, Descriptor>().generate(unserializer));
        serializerPos = unserializer.getCurrentPos();
    }
}

/* *************** class RestParticle2D ************************************ */

template <typename T, template <typename U> class Descriptor>
int RestParticle2D<T, Descriptor>::id =
    meta::registerGenericParticle2D<T, Descriptor, RestParticle2D<T, Descriptor> >("Rest");

template <typename T, template <typename U> class Descriptor>
RestParticle2D<T, Descriptor>::RestParticle2D()
{ }

template <typename T, template <typename U> class Descriptor>
RestParticle2D<T, Descriptor>::RestParticle2D(plint tag_, Array<T, 2> const &position) :
    Particle2D<T, Descriptor>(tag_, position)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle2D<T, Descriptor>::velocityToParticle(
    TensorField2D<T, 2> &velocityField, T scaling)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle2D<T, Descriptor>::velocityToParticle(NTensorField2D<T> &velocityField, T scaling)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle2D<T, Descriptor>::rhoBarJtoParticle(
    NTensorField2D<T> &rhoBarJfield, bool velIsJ, T scaling)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle2D<T, Descriptor>::fluidToParticle(BlockLattice2D<T, Descriptor> &fluid, T scaling)
{ }

template <typename T, template <typename U> class Descriptor>
void RestParticle2D<T, Descriptor>::advance()
{ }

template <typename T, template <typename U> class Descriptor>
int RestParticle2D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
RestParticle2D<T, Descriptor> *RestParticle2D<T, Descriptor>::clone() const
{
    return new RestParticle2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
bool RestParticle2D<T, Descriptor>::getVector(plint whichVector, Array<T, 2> &vector) const
{
    return Particle2D<T, Descriptor>::getVector(whichVector, vector);
}

/* *************** class VerletParticle2D ************************************ */

template <typename T, template <typename U> class Descriptor>
int VerletParticle2D<T, Descriptor>::id =
    meta::registerGenericParticle2D<T, Descriptor, VerletParticle2D<T, Descriptor> >("Verlet");

template <typename T, template <typename U> class Descriptor>
VerletParticle2D<T, Descriptor>::VerletParticle2D() :
    v(T(), T()), vHalfTime(T(), T()), a(T(), T()), fluidCompliance((T)1.), rho((T)1.), invRho((T)1.)
{ }

template <typename T, template <typename U> class Descriptor>
VerletParticle2D<T, Descriptor>::VerletParticle2D(plint tag_, Array<T, 2> const &position) :
    Particle2D<T, Descriptor>(tag_, position),
    v(T(), T()),
    vHalfTime(T(), T()),
    a(T(), T()),
    fluidCompliance((T)1.),
    rho((T)1.),
    invRho((T)1.)
{ }

template <typename T, template <typename U> class Descriptor>
void VerletParticle2D<T, Descriptor>::velocityToParticle(
    TensorField2D<T, 2> &velocityField, T scaling)
{
    Array<T, 2> position(this->getPosition());
    std::vector<Dot2D> pos(4);
    std::vector<T> weights(4);
    linearInterpolationCoefficients(velocityField, position, pos, weights);
    Array<T, 2> fluidVelocity;
    fluidVelocity.resetToZero();
    for (plint iCell = 0; iCell < 4; ++iCell) {
        fluidVelocity += weights[iCell] * velocityField.get(pos[iCell].x, pos[iCell].y) * scaling;
    }

    Array<T, 2> force((fluidVelocity - this->get_v()) * fluidCompliance);
    this->set_a(force / this->get_rho());
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle2D<T, Descriptor>::velocityToParticle(
    NTensorField2D<T> &velocityField, T scaling)
{
    Array<T, 2> position(this->getPosition());
    std::vector<Dot2D> pos(4);
    std::vector<T> weights(4);
    linearInterpolationCoefficients(velocityField, position, pos, weights);
    Array<T, 2> fluidVelocity;
    fluidVelocity.resetToZero();
    for (plint iCell = 0; iCell < 4; ++iCell) {
        T *data = velocityField.get(pos[iCell].x, pos[iCell].y);
        fluidVelocity[0] += weights[iCell] * scaling * data[0];
        fluidVelocity[1] += weights[iCell] * scaling * data[1];
    }

    Array<T, 2> force((fluidVelocity - this->get_v()) * fluidCompliance);
    this->set_a(force / this->get_rho());
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle2D<T, Descriptor>::rhoBarJtoParticle(
    NTensorField2D<T> &rhoBarJfield, bool velIsJ, T scaling)
{
    Array<T, 2> position(this->getPosition());
    std::vector<Dot2D> pos(4);
    std::vector<T> weights(4);
    linearInterpolationCoefficients(rhoBarJfield, position, pos, weights);
    Array<T, 2> j;
    j.resetToZero();
    T rhoBar = T();
    for (plint iCell = 0; iCell < 4; ++iCell) {
        T *data = rhoBarJfield.get(pos[iCell].x, pos[iCell].y);
        j.add_from_cArray(data + 1, weights[iCell]);
        rhoBar += weights[iCell] * (*data);
    }
    Array<T, 2> fluidVelocity;
    if (velIsJ) {
        fluidVelocity = j * scaling;
    } else {
        fluidVelocity = j * scaling * Descriptor<T>::invRho(rhoBar);
    }

    Array<T, 2> force((fluidVelocity - this->get_v()) * fluidCompliance);
    this->set_a(force / this->get_rho());
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle2D<T, Descriptor>::fluidToParticle(
    BlockLattice2D<T, Descriptor> &fluid, T scaling)
{
    Array<T, 2> position(this->getPosition());
    std::vector<Dot2D> pos(4);
    std::vector<T> weights(4);
    linearInterpolationCoefficients(fluid, position, pos, weights);
    Array<T, 2> tmpVel;
    Array<T, 2> fluidVelocity;
    fluidVelocity.resetToZero();
    for (plint iCell = 0; iCell < 4; ++iCell) {
        if (contained(pos[iCell].x, pos[iCell].y, fluid.getBoundingBox())) {
            fluid.get(pos[iCell].x, pos[iCell].y).computeVelocity(tmpVel);
        } else {
            tmpVel.resetToZero();
        }
        fluidVelocity += weights[iCell] * tmpVel * scaling;
    }

    Array<T, 2> force((fluidVelocity - this->get_v()) * fluidCompliance);
    this->set_a(force / this->get_rho());
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle2D<T, Descriptor>::advance()
{
    vHalfTime = v + (T)0.5 * a;
    this->getPosition() += vHalfTime;
}

template <typename T, template <typename U> class Descriptor>
int VerletParticle2D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle2D<T, Descriptor>::reset(Array<T, 2> const &position_)
{
    Particle2D<T, Descriptor>::reset(position_);
    v.resetToZero();
    vHalfTime.resetToZero();
    a.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle2D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    Particle2D<T, Descriptor>::serialize(serializer);
    serializer.addValues<T, 2>(v);
    serializer.addValues<T, 2>(vHalfTime);
    serializer.addValues<T, 2>(a);
    serializer.addValue<T>(fluidCompliance);
    serializer.addValue<T>(rho);
    serializer.addValue<T>(invRho);
}

template <typename T, template <typename U> class Descriptor>
void VerletParticle2D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    Particle2D<T, Descriptor>::unserialize(unserializer);
    unserializer.readValues<T, 2>(v);
    unserializer.readValues<T, 2>(vHalfTime);
    unserializer.readValues<T, 2>(a);
    unserializer.readValue<T>(fluidCompliance);
    unserializer.readValue<T>(rho);
    unserializer.readValue<T>(invRho);
}

template <typename T, template <typename U> class Descriptor>
VerletParticle2D<T, Descriptor> *VerletParticle2D<T, Descriptor>::clone() const
{
    return new VerletParticle2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
bool VerletParticle2D<T, Descriptor>::getVector(plint whichVector, Array<T, 2> &vector) const
{
    if (whichVector == 0) {
        vector = get_v();
        return true;
    } else if (whichVector == 1) {
        vector = get_a();
        return true;
    }
    return Particle2D<T, Descriptor>::getVector(whichVector, vector);
}

template <typename T, template <typename U> class Descriptor>
bool VerletParticle2D<T, Descriptor>::getScalar(plint whichScalar, T &scalar) const
{
    if (whichScalar == 0) {
        scalar = rho;
        return true;
    } else if (whichScalar == 1) {
        scalar = fluidCompliance;
        return true;
    }
    return Particle2D<T, Descriptor>::getScalar(whichScalar, scalar);
}

/* *************** class VerletGeryaParticle2D ************************************ */
template <typename T, template <typename U> class Descriptor>
int VerletGeryaParticle2D<T, Descriptor>::id =
    meta::registerGenericParticle2D<T, Descriptor, VerletGeryaParticle2D<T, Descriptor> >(
        "VerletGerya");

template <typename T, template <typename U> class Descriptor>
VerletGeryaParticle2D<T, Descriptor>::VerletGeryaParticle2D() :
    VerletParticle2D<T, Descriptor>(), fluidConcentration((T)0)
{ }

template <typename T, template <typename U> class Descriptor>
VerletGeryaParticle2D<T, Descriptor>::VerletGeryaParticle2D(
    plint tag_, Array<T, 2> const &position_, T fluidConcentration_) :
    VerletParticle2D<T, Descriptor>(tag_, position_), fluidConcentration(fluidConcentration_)
{ }

template <typename T, template <typename U> class Descriptor>
int VerletGeryaParticle2D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
VerletGeryaParticle2D<T, Descriptor> *VerletGeryaParticle2D<T, Descriptor>::clone() const
{
    return new VerletGeryaParticle2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void VerletGeryaParticle2D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    VerletParticle2D<T, Descriptor>::serialize(serializer);
    serializer.addValue<T>(fluidConcentration);
}

template <typename T, template <typename U> class Descriptor>
void VerletGeryaParticle2D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    VerletParticle2D<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue<T>(fluidConcentration);
}

template <typename T, template <typename U> class Descriptor>
void VerletGeryaParticle2D<T, Descriptor>::reset(Array<T, 2> const &position_)
{
    VerletParticle2D<T, Descriptor>::reset(position_);
    fluidConcentration = T();
}

template <typename T, template <typename U> class Descriptor>
T VerletGeryaParticle2D<T, Descriptor>::getFluidConcentration()
{
    return fluidConcentration;
}

template <typename T, template <typename U> class Descriptor>
void VerletGeryaParticle2D<T, Descriptor>::setFluidConcentration(T fluidConcentration_)
{
    fluidConcentration = fluidConcentration_;
}

}  // namespace plb

#endif  // PARTICLE_2D_H
