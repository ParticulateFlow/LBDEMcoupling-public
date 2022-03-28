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

#ifndef VISUAL_PARTICLE_2D_HH
#define VISUAL_PARTICLE_2D_HH

#include "core/globalDefs.h"
#include "particles/visualParticle2D.h"

namespace plb {

/* *************** class VisualParticle2D ************************************ */

template <typename T, template <typename U> class Descriptor>
int VisualParticle2D<T, Descriptor>::id =
    meta::registerVisualParticle2D<T, Descriptor, VisualParticle2D<T, Descriptor> >("Visual");

template <typename T, template <typename U> class Descriptor>
VisualParticle2D<T, Descriptor>::VisualParticle2D()
{ }

template <typename T, template <typename U> class Descriptor>
VisualParticle2D<T, Descriptor>::VisualParticle2D(plint tag_, Array<T, 2> const &position) :
    Particle2D<T, Descriptor>(tag_, position)
{ }

template <typename T, template <typename U> class Descriptor>
VisualParticle2D<T, Descriptor>::VisualParticle2D(
    plint tag_, Array<T, 2> const &position, std::vector<T> const &scalars_,
    std::vector<Array<T, 2> > const &vectors_) :
    Particle2D<T, Descriptor>(tag_, position), scalars(scalars_), vectors(vectors_)
{ }

template <typename T, template <typename U> class Descriptor>
int VisualParticle2D<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void VisualParticle2D<T, Descriptor>::reset(Array<T, 2> const &position_)
{
    Particle2D<T, Descriptor>::reset(position_);
    std::fill(scalars.begin(), scalars.end(), T());
    std::fill(vectors.begin(), vectors.end(), Array<T, 2>(T(), T(), T()));
}

template <typename T, template <typename U> class Descriptor>
void VisualParticle2D<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    Particle2D<T, Descriptor>::serialize(serializer);
    pluint scalarSize = scalars.size();
    serializer.addValue(scalarSize);
    for (pluint iScalar = 0; iScalar < scalars.size(); ++iScalar) {
        serializer.addValue(scalars[iScalar]);
    }
    pluint vectorSize = vectors.size();
    serializer.addValue(vectorSize);
    for (pluint iVector = 0; iVector < vectors.size(); ++iVector) {
        serializer.addValues<T, 2>(vectors[iVector]);
    }
}

template <typename T, template <typename U> class Descriptor>
void VisualParticle2D<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    Particle2D<T, Descriptor>::unserialize(unserializer);
    pluint scalarSize;
    unserializer.readValue(scalarSize);
    PLB_ASSERT(scalarSize == (pluint)scalars.size());
    for (pluint iScalar = 0; iScalar < scalars.size(); ++iScalar) {
        unserializer.readValue(scalars[iScalar]);
    }
    pluint vectorSize;
    unserializer.readValue(vectorSize);
    PLB_ASSERT(vectorSize == (pluint)vectors.size());
    for (pluint iVector = 0; iVector < vectors.size(); ++iVector) {
        unserializer.readValues<T, 2>(vectors[iVector]);
    }
}

template <typename T, template <typename U> class Descriptor>
VisualParticle2D<T, Descriptor> *VisualParticle2D<T, Descriptor>::clone() const
{
    return new VisualParticle2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
bool VisualParticle2D<T, Descriptor>::getVector(plint whichVector, Array<T, 2> &vector) const
{
    if (whichVector < (plint)vectors.size()) {
        vector = vectors[whichVector];
        return true;
    }
    return Particle2D<T, Descriptor>::getVector(whichVector, vector);
}

template <typename T, template <typename U> class Descriptor>
bool VisualParticle2D<T, Descriptor>::getScalar(plint whichScalar, T &scalar) const
{
    if (whichScalar < (plint)scalars.size()) {
        scalar = scalars[whichScalar];
        return true;
    }
    return Particle2D<T, Descriptor>::getScalar(whichScalar, scalar);
}

template <typename T, template <typename U> class Descriptor>
bool VisualParticle2D<T, Descriptor>::setScalar(plint whichScalar, T scalar)
{
    if (whichScalar < (plint)scalars.size()) {
        scalars[whichScalar] = scalar;
        return true;
    } else {
        return false;
    }
}

template <typename T, template <typename U> class Descriptor>
bool VisualParticle2D<T, Descriptor>::setScalars(std::vector<T> const &scalars_)
{
    scalars = scalars_;
    return true;
}

template <typename T, template <typename U> class Descriptor>
bool VisualParticle2D<T, Descriptor>::setVectors(std::vector<Array<T, 2> > const &vectors_)
{
    vectors = vectors_;
    return true;
}

template <typename T, template <typename U> class Descriptor>
bool VisualParticle2D<T, Descriptor>::setVector(plint whichVector, Array<T, 2> const &vector)
{
    if (whichVector < (plint)vectors.size()) {
        vectors[whichVector] = vector;
        return true;
    } else {
        return false;
    }
}

// /* *************** class MappingParticle2D ************************************ */
//
// template<typename T, template<typename U> class Descriptor>
// MappingParticle2D<T,Descriptor>::MappingParticle2D()
//     : PointParticle2D<T,Descriptor>(),
//       stage(0),
//       pos1(T(),T(),T()),
//       pos2(T(),T(),T())
// { }
//
// template<typename T, template<typename U> class Descriptor>
// MappingParticle2D<T,Descriptor>::MappingParticle2D(plint tag_, Array<T,2> const& position_,
// Array<T,2> const& velocity_)
//     : PointParticle2D<T,Descriptor>(tag_, position_, velocity_),
//       stage(0),
//       pos1(T(),T(),T()),
//       pos2(T(),T(),T())
// { }
//
// template<typename T, template<typename U> class Descriptor>
// void MappingParticle2D<T,Descriptor>::advance() {
//     // Dynamics of a mapping particle.
//     // Stage 0: The particle is behind the first surface and is flowing forward.
//     // Stage 1: The particle is behind the second surface and is flowing forward.
//     // Stage 2: The particle has crossed the second surface and is
//     //          now walking back along a straight line to the first surface.
//     // Stage 3: The particle is back at the first surface and does nothing any more.
//
//     // Forward flowing.
//     if (stage==0 || stage==1) {
//         this->getPosition() += this->getVelocity();
//     }
//     // Backward walking.
//     else if (stage==2) {
//         Array<T,2> homeArrow = pos1-this->getPosition();
//         T l = norm(homeArrow);
//         T step=0.5;
//         if (l<step) {
//             stage=3;
//             this->getPosition() = pos1;
//         }
//         else {
//             this->getPosition() += homeArrow *step/l;
//         }
//     }
//     // Crossed first surface?
//     if (stage==0 && crossedSurface1()) {
//         stage=1;
//         pos1 = this->getPosition();
//     }
//     // Crossed second surface?
//     if (stage==1 && crossedSurface2()) {
//         stage=2;
//         pos2 = this->getPosition();
//     }
// }
//
// template<typename T, template<typename U> class Descriptor>
// void MappingParticle2D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
// {
//     PointParticle2D<T,Descriptor>::serialize(serializer);
//     serializer.addValue(stage);
//     serializer.addValues<T,2>(pos1);
//     serializer.addValues<T,2>(pos2);
// }
//
// template<typename T, template<typename U> class Descriptor>
// void MappingParticle2D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
// {
//     PointParticle2D<T,Descriptor>::unserialize(unserializer);
//     unserializer.readValue(stage);
//     unserializer.readValues<T,2>(pos1);
//     unserializer.readValues<T,2>(pos2);
// }
//
// template<typename T, template<typename U> class Descriptor>
// void MappingParticle2D<T,Descriptor>::rescale(int dxScale, int dtScale) {
//     PointParticle2D<T,Descriptor>::rescale(dxScale, dtScale);
//     T posDimDx = 1.;
//     T posDimDt = 0.;
//     T posSaleFactor = scaleFromReference(dxScale, posDimDx, dtScale, posDimDt);
//     pos1 *= posSaleFactor;
//     pos2 *= posSaleFactor;
// }
//
//
// /* *************** class MappingParticleZslice2D ************************************ */
//
// template<typename T, template<typename U> class Descriptor>
// int MappingParticleZslice2D<T,Descriptor>::id =
//         meta::registerGenericParticle2D<T,Descriptor,MappingParticleZslice2D<T,Descriptor>
//         >("MappingZslice");
//
// template<typename T, template<typename U> class Descriptor>
// MappingParticleZslice2D<T,Descriptor>::MappingParticleZslice2D()
//     : MappingParticle2D<T,Descriptor>(),
//       zSlice1(T()), zSlice2(T())
// { }
//
// template<typename T, template<typename U> class Descriptor>
// MappingParticleZslice2D<T,Descriptor>::MappingParticleZslice2D (
//         plint tag_, Array<T,2> const& position_, Array<T,2> const& velocity_, T zSlice1_, T
//         zSlice2_)
//     : MappingParticle2D<T,Descriptor>(tag_, position_, velocity_),
//       zSlice1(zSlice1_),
//       zSlice2(zSlice2_)
// { }
//
// template<typename T, template<typename U> class Descriptor>
// int MappingParticleZslice2D<T,Descriptor>::getId() const
// {
//     return id;
// }
//
// template<typename T, template<typename U> class Descriptor>
// bool MappingParticleZslice2D<T,Descriptor>::crossedSurface1() const {
//     if (zSlice2>zSlice1) {
//         return this->getPosition()[2] > zSlice1;
//     }
//     else {
//         return this->getPosition()[2] < zSlice1;
//     }
// }
//
// template<typename T, template<typename U> class Descriptor>
// bool MappingParticleZslice2D<T,Descriptor>::crossedSurface2() const {
//     if (zSlice2>zSlice1) {
//         return this->getPosition()[2] > zSlice2;
//     }
//     else {
//         return this->getPosition()[2] < zSlice2;
//     }
// }
//
// template<typename T, template<typename U> class Descriptor>
// T MappingParticleZslice2D<T,Descriptor>::getSurfaceDistance() const {
//     return std::fabs(zSlice2 - zSlice1);
// }
//
// template<typename T, template<typename U> class Descriptor>
// MappingParticleZslice2D<T,Descriptor>* MappingParticleZslice2D<T,Descriptor>::clone() const {
//     return new MappingParticleZslice2D<T,Descriptor>(*this);
// }
//
//
// template<typename T, template<typename U> class Descriptor>
// void MappingParticleZslice2D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
// {
//     MappingParticle2D<T,Descriptor>::serialize(serializer);
//     serializer.addValue(zSlice1);
//     serializer.addValue(zSlice2);
// }
//
// template<typename T, template<typename U> class Descriptor>
// void MappingParticleZslice2D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
// {
//     MappingParticle2D<T,Descriptor>::unserialize(unserializer);
//     unserializer.readValue(zSlice1);
//     unserializer.readValue(zSlice2);
// }
//
//
// template<typename T, template<typename U> class Descriptor>
// void MappingParticleZslice2D<T,Descriptor>::rescale(int dxScale, int dtScale) {
//     MappingParticle2D<T,Descriptor>::rescale(dxScale, dtScale);
//     T posDimDx = 1.;
//     T posDimDt = 0.;
//     T posSaleFactor = scaleFromReference(dxScale, posDimDx, dtScale, posDimDt);
//     zSlice1 *= posSaleFactor;
//     zSlice2 *= posSaleFactor;
// }

}  // namespace plb

#endif  // VISUAL_PARTICLE_2D_H
