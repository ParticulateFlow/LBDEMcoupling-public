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

#ifndef OFF_LATTICE_BOUNDARY_PROFILES_3D_HH
#define OFF_LATTICE_BOUNDARY_PROFILES_3D_HH

#include "core/functions.h"
#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/offLatticeBoundaryProfiles3D.h"

namespace plb {

template <typename T>
struct DefaultWallProfile3D<T, Array<T, 3> > {
    BoundaryProfile3D<T, Array<T, 3> > *generate()
    {
        return new NoSlipProfile3D<T>();
    }
};

template <typename T>
struct DefaultWallProfile3D<T, Array<T, 2> > {
    BoundaryProfile3D<T, Array<T, 2> > *generate()
    {
        return new ScalarNeumannProfile3D<T>();
    }
};

/********** NoSlipProfile3D ********************************************/

template <typename T>
void NoSlipProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{ }

template <typename T>
void NoSlipProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void NoSlipProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    data.resetToZero();
    bdType = OffBoundary::dirichlet;
}

template <typename T>
NoSlipProfile3D<T> *NoSlipProfile3D<T>::clone() const
{
    return new NoSlipProfile3D<T>(*this);
}

/********** FreeSlipProfile3D ********************************************/

template <typename T>
void FreeSlipProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{ }

template <typename T>
void FreeSlipProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void FreeSlipProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    data.resetToZero();
    bdType = OffBoundary::freeSlip;
}

template <typename T>
FreeSlipProfile3D<T> *FreeSlipProfile3D<T>::clone() const
{
    return new FreeSlipProfile3D<T>(*this);
}

/********** PoiseuilleProfile3D ********************************************/

template <typename T>
PoiseuilleProfile3D<T>::PoiseuilleProfile3D(T uAverage_) :
    uAverage(uAverage_), normal(T(), T(), T()), center(T(), T(), T()), radius((T)1)
{ }

template <typename T>
void PoiseuilleProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{
    normal = normal_;
}

template <typename T>
void PoiseuilleProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{
    center = center_;
    radius = radius_;
}

template <typename T>
void PoiseuilleProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::dirichlet;
    Array<T, 3> radial = pos - center;
    T r = norm(radial) / radius;
    if (r <= (T)1.) {
        data = 2 * uAverage * (1 - r * r) * (-normal);
    } else {
        data.resetToZero();
    }
}

template <typename T>
PoiseuilleProfile3D<T> *PoiseuilleProfile3D<T>::clone() const
{
    return new PoiseuilleProfile3D<T>(*this);
}

/********** IncreasingPoiseuilleProfile3D ********************************************/

template <typename T, template <typename U> class Descriptor>
IncreasingPoiseuilleProfile3D<T, Descriptor>::IncreasingPoiseuilleProfile3D(
    T uAverage_, plint maxT_) :
    uAverage(uAverage_), maxT(maxT_), normal(T(), T(), T()), center(T(), T(), T()), radius((T)1)
{ }

template <typename T, template <typename U> class Descriptor>
void IncreasingPoiseuilleProfile3D<T, Descriptor>::setNormal(Array<T, 3> const &normal_)
{
    normal = normal_;
}

template <typename T, template <typename U> class Descriptor>
void IncreasingPoiseuilleProfile3D<T, Descriptor>::defineCircularShape(
    Array<T, 3> const &center_, T radius_)
{
    center = center_;
    radius = radius_;
}

template <typename T, template <typename U> class Descriptor>
void IncreasingPoiseuilleProfile3D<T, Descriptor>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::dirichlet;
    BlockLattice3D<T, Descriptor> const *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> const *>(argument);
    plint t = lattice->getTimeCounter().getTime();
    T signal = util::sinIncreasingFunction<T>((T)t, (T)maxT);
    Array<T, 3> radial = pos - center;
    T r = norm(radial) / radius;
    if (r <= (T)1.) {
        data = signal * 2 * uAverage * (1 - r * r) * (-normal);
    } else {
        data.resetToZero();
    }
}

template <typename T, template <typename U> class Descriptor>
IncreasingPoiseuilleProfile3D<T, Descriptor> *IncreasingPoiseuilleProfile3D<T, Descriptor>::clone()
    const
{
    return new IncreasingPoiseuilleProfile3D<T, Descriptor>(*this);
}

/********** IncreasingVelocityProfile3D ********************************************/

template <typename T, template <typename U> class Descriptor>
IncreasingVelocityProfile3D<T, Descriptor>::IncreasingVelocityProfile3D(
    Array<T, 3> const &u_, plint maxT_) :
    u(u_), maxT(maxT_)
{ }

template <typename T, template <typename U> class Descriptor>
void IncreasingVelocityProfile3D<T, Descriptor>::setNormal(Array<T, 3> const &normal_)
{ }

template <typename T, template <typename U> class Descriptor>
void IncreasingVelocityProfile3D<T, Descriptor>::defineCircularShape(
    Array<T, 3> const &center_, T radius_)
{ }

template <typename T, template <typename U> class Descriptor>
void IncreasingVelocityProfile3D<T, Descriptor>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::dirichlet;
    BlockLattice3D<T, Descriptor> const *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> const *>(argument);
    plint t = lattice->getTimeCounter().getTime();
    T signal = util::sinIncreasingFunction<T>((T)t, (T)maxT);
    data = u * signal;
}

template <typename T, template <typename U> class Descriptor>
IncreasingVelocityProfile3D<T, Descriptor> *IncreasingVelocityProfile3D<T, Descriptor>::clone()
    const
{
    return new IncreasingVelocityProfile3D<T, Descriptor>(*this);
}

/********** IncreasingVelocityPlugProfile3D ********************************************/

template <typename T, template <typename U> class Descriptor>
IncreasingVelocityPlugProfile3D<T, Descriptor>::IncreasingVelocityPlugProfile3D(
    T uMax_, plint maxT_) :
    uMax(uMax_), maxT(maxT_), normal(T(), T(), T())
{ }

template <typename T, template <typename U> class Descriptor>
void IncreasingVelocityPlugProfile3D<T, Descriptor>::setNormal(Array<T, 3> const &normal_)
{
    normal = normal_;
}

template <typename T, template <typename U> class Descriptor>
void IncreasingVelocityPlugProfile3D<T, Descriptor>::defineCircularShape(
    Array<T, 3> const &center_, T radius_)
{ }

template <typename T, template <typename U> class Descriptor>
void IncreasingVelocityPlugProfile3D<T, Descriptor>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::dirichlet;
    BlockLattice3D<T, Descriptor> const *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> const *>(argument);
    plint t = lattice->getTimeCounter().getTime();
    T signal = util::sinIncreasingFunction<T>((T)t, (T)maxT);
    data = -signal * uMax * normal;
}

template <typename T, template <typename U> class Descriptor>
IncreasingVelocityPlugProfile3D<T, Descriptor>
    *IncreasingVelocityPlugProfile3D<T, Descriptor>::clone() const
{
    return new IncreasingVelocityPlugProfile3D<T, Descriptor>(*this);
}

/********** TimeDependentVelocityProfile3D ********************************************/

template <typename T, template <typename U> class Descriptor>
TimeDependentVelocityProfile3D<T, Descriptor>::TimeDependentVelocityProfile3D(
    util::TimeDependentFunction<T, 3> *velocity_) :
    velocity(velocity_)
{ }

template <typename T, template <typename U> class Descriptor>
TimeDependentVelocityProfile3D<T, Descriptor>::TimeDependentVelocityProfile3D(
    TimeDependentVelocityProfile3D<T, Descriptor> const &rhs) :
    velocity(rhs.velocity->clone())
{ }

template <typename T, template <typename U> class Descriptor>
TimeDependentVelocityProfile3D<T, Descriptor>
    &TimeDependentVelocityProfile3D<T, Descriptor>::operator=(
        TimeDependentVelocityProfile3D<T, Descriptor> const &rhs)
{
    TimeDependentVelocityProfile3D<T, Descriptor>(rhs).swap(*this);
}

template <typename T, template <typename U> class Descriptor>
void TimeDependentVelocityProfile3D<T, Descriptor>::swap(
    TimeDependentVelocityProfile3D<T, Descriptor> &rhs)
{
    std::swap(velocity, rhs.velocity);
}

template <typename T, template <typename U> class Descriptor>
void TimeDependentVelocityProfile3D<T, Descriptor>::setNormal(Array<T, 3> const &normal_)
{ }

template <typename T, template <typename U> class Descriptor>
void TimeDependentVelocityProfile3D<T, Descriptor>::defineCircularShape(
    Array<T, 3> const &center_, T radius_)
{ }

template <typename T, template <typename U> class Descriptor>
void TimeDependentVelocityProfile3D<T, Descriptor>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::dirichlet;
    BlockLattice3D<T, Descriptor> const *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> const *>(argument);
    plint t = lattice->getTimeCounter().getTime();
    data = (*velocity)(t);
}

template <typename T, template <typename U> class Descriptor>
TimeDependentVelocityProfile3D<T, Descriptor>
    *TimeDependentVelocityProfile3D<T, Descriptor>::clone() const
{
    return new TimeDependentVelocityProfile3D<T, Descriptor>(*this);
}

/********** OscillatingPoiseuilleProfile3D ********************************************/

template <typename T, template <typename U> class Descriptor>
OscillatingPoiseuilleProfile3D<T, Descriptor>::OscillatingPoiseuilleProfile3D(
    T minUave_, T maxUave_, T period_) :
    minUave(minUave_),
    maxUave(maxUave_),
    period(period_),
    normal(T(), T(), T()),
    center(T(), T(), T()),
    radius((T)1)
{ }

template <typename T, template <typename U> class Descriptor>
void OscillatingPoiseuilleProfile3D<T, Descriptor>::setNormal(Array<T, 3> const &normal_)
{
    normal = normal_;
}

template <typename T, template <typename U> class Descriptor>
void OscillatingPoiseuilleProfile3D<T, Descriptor>::defineCircularShape(
    Array<T, 3> const &center_, T radius_)
{
    center = center_;
    radius = radius_;
}

template <typename T, template <typename U> class Descriptor>
void OscillatingPoiseuilleProfile3D<T, Descriptor>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::dirichlet;
    static const T pi = (T)4. * std::atan((T)1.);
    BlockLattice3D<T, Descriptor> const *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> const *>(argument);
    PLB_ASSERT(lattice);
    T t = (T)lattice->getTimeCounter().getTime();
    T signal = (std::sin((T)2. * pi * t / period) + (T)1.) * (T)0.5;

    Array<T, 3> radial = pos - center;
    T r = norm(radial) / radius;
    if (r <= (T)1.) {
        data = ((T)2 * (minUave + (maxUave - minUave) * signal)) * ((T)1 - r * r) * (-normal);
    } else {
        data.resetToZero();
    }
}

template <typename T, template <typename U> class Descriptor>
OscillatingPoiseuilleProfile3D<T, Descriptor>
    *OscillatingPoiseuilleProfile3D<T, Descriptor>::clone() const
{
    return new OscillatingPoiseuilleProfile3D<T, Descriptor>(*this);
}

/********** ConstantVelocityProfile3D ********************************************/

template <typename T>
ConstantVelocityProfile3D<T>::ConstantVelocityProfile3D(Array<T, 3> const &u_) : u(u_)
{ }

template <typename T>
void ConstantVelocityProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{ }

template <typename T>
void ConstantVelocityProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void ConstantVelocityProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::dirichlet;
    data = u;
}

template <typename T>
ConstantVelocityProfile3D<T> *ConstantVelocityProfile3D<T>::clone() const
{
    return new ConstantVelocityProfile3D<T>(*this);
}

/********** VelocityPlugProfile3D ********************************************/

template <typename T>
VelocityPlugProfile3D<T>::VelocityPlugProfile3D(T uMax_) : uMax(uMax_), normal(T(), T(), T())
{ }

template <typename T>
void VelocityPlugProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{
    normal = normal_;
}

template <typename T>
void VelocityPlugProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void VelocityPlugProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::dirichlet;
    data = -uMax * normal;
}

template <typename T>
VelocityPlugProfile3D<T> *VelocityPlugProfile3D<T>::clone() const
{
    return new VelocityPlugProfile3D<T>(*this);
}

/********** ExactRotationalVelocityProfile3D ********************************************/

template <typename T>
ExactRotationalVelocityProfile3D<T>::ExactRotationalVelocityProfile3D(
    Array<T, 3> const &angularVelocity_, Array<T, 3> const &pointOnRotationAxis_) :
    angularVelocity(angularVelocity_), pointOnRotationAxis(pointOnRotationAxis_)
{ }

template <typename T>
void ExactRotationalVelocityProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{ }

template <typename T>
void ExactRotationalVelocityProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void ExactRotationalVelocityProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::dirichlet;
    data = getExactRotationalVelocity(pos, angularVelocity, pointOnRotationAxis);
}

template <typename T>
ExactRotationalVelocityProfile3D<T> *ExactRotationalVelocityProfile3D<T>::clone() const
{
    return new ExactRotationalVelocityProfile3D<T>(*this);
}

/********** NeumannBoundaryProfile3D ******************************************/

template <typename T>
void NeumannBoundaryProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{
    normal = normal_;
}

template <typename T>
void NeumannBoundaryProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void NeumannBoundaryProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::neumann;
    // In a Neumann condition, the velocity will need to be
    //   computed in the actual off-lattice boundary algorithm.
    //   Nothing to be done here.
    data.resetToZero();
}

template <typename T>
NeumannBoundaryProfile3D<T> *NeumannBoundaryProfile3D<T>::clone() const
{
    return new NeumannBoundaryProfile3D<T>(*this);
}

/********** DensityNeumannBoundaryProfile3D ******************************************/

template <typename T>
DensityNeumannBoundaryProfile3D<T>::DensityNeumannBoundaryProfile3D(T rho_) : rho(rho_)
{ }

template <typename T>
void DensityNeumannBoundaryProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{ }

template <typename T>
void DensityNeumannBoundaryProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void DensityNeumannBoundaryProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 3> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::densityNeumann;
    // In a Neumann condition, the velocity will need to be
    //   computed in the actual off-lattice boundary algorithm.
    //   Nothing to be done here.
    data.resetToZero();
    data[0] = rho;
}

template <typename T>
DensityNeumannBoundaryProfile3D<T> *DensityNeumannBoundaryProfile3D<T>::clone() const
{
    return new DensityNeumannBoundaryProfile3D<T>(*this);
}

/********** ScalarDirichletProfile3D ******************************************/

template <typename T>
ScalarDirichletProfile3D<T>::ScalarDirichletProfile3D(T value_) : value(value_)
{ }

template <typename T>
void ScalarDirichletProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{ }

template <typename T>
void ScalarDirichletProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void ScalarDirichletProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 2> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::dirichlet;
    data[0] = value;
    data[1] = T();  // Second argument is not used here.
}

template <typename T>
ScalarDirichletProfile3D<T> *ScalarDirichletProfile3D<T>::clone() const
{
    return new ScalarDirichletProfile3D<T>(*this);
}

/********** ScalarNeumannProfile3D ******************************************/

template <typename T>
void ScalarNeumannProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{
    normal = normal_;
}

template <typename T>
void ScalarNeumannProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void ScalarNeumannProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 2> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::neumann;
    // Neumann means zero-gradient, the two arguments are unused.
    data[0] = T();
    data[1] = T();
}

template <typename T>
ScalarNeumannProfile3D<T> *ScalarNeumannProfile3D<T>::clone() const
{
    return new ScalarNeumannProfile3D<T>(*this);
}

/********** ScalarFluxProfile3D ******************************************/

template <typename T>
ScalarFluxProfile3D<T>::ScalarFluxProfile3D(T gradVal_) : gradVal(gradVal_)
{ }

template <typename T>
void ScalarFluxProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{
    normal = normal_;
}

template <typename T>
void ScalarFluxProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void ScalarFluxProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 2> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::flux;
    data[0] = gradVal;
    data[1] = T();
}

template <typename T>
ScalarFluxProfile3D<T> *ScalarFluxProfile3D<T>::clone() const
{
    return new ScalarFluxProfile3D<T>(*this);
}

/********** ScalarIsolationProfile3D ******************************************/

template <typename T>
ScalarIsolationProfile3D<T>::ScalarIsolationProfile3D(T asymptoticRho_, T kappa_) :
    asymptoticRho(asymptoticRho_), kappa(kappa_)
{ }

template <typename T>
void ScalarIsolationProfile3D<T>::setNormal(Array<T, 3> const &normal_)
{
    normal = normal_;
}

template <typename T>
void ScalarIsolationProfile3D<T>::defineCircularShape(Array<T, 3> const &center_, T radius_)
{ }

template <typename T>
void ScalarIsolationProfile3D<T>::getData(
    Array<T, 3> const &pos, plint id, AtomicBlock3D const *argument, Array<T, 2> &data,
    OffBoundary::Type &bdType) const
{
    bdType = OffBoundary::isolation;
    data[0] = asymptoticRho;
    data[1] = kappa;
}

template <typename T>
ScalarIsolationProfile3D<T> *ScalarIsolationProfile3D<T>::clone() const
{
    return new ScalarIsolationProfile3D<T>(*this);
}

}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_PROFILES_3D_HH
