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

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef BOUNCE_BACK_MODELS_HH
#define BOUNCE_BACK_MODELS_HH

#include <algorithm>
#include <limits>

#include "boundaryCondition/bounceBackModels.h"
#include "core/array.h"
#include "core/cell.h"
#include "core/latticeStatistics.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

/* *************** Class MomentumExchangeBounceBack ************************* */

template <typename T, template <typename U> class Descriptor>
MomentumExchangeBounceBack<T, Descriptor>::MomentumExchangeBounceBack(
    Array<plint, Descriptor<T>::d> forceIds_, T rho_) :
    fluidDirections(), forceIds(forceIds_), rho(rho_)
{ }

template <typename T, template <typename U> class Descriptor>
MomentumExchangeBounceBack<T, Descriptor>::MomentumExchangeBounceBack(
    HierarchicUnserializer &unserializer)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
MomentumExchangeBounceBack<T, Descriptor> *MomentumExchangeBounceBack<T, Descriptor>::clone() const
{
    return new MomentumExchangeBounceBack<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int MomentumExchangeBounceBack<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, MomentumExchangeBounceBack<T, Descriptor> >(
        "MomentumExchangeBounceBack");

template <typename T, template <typename U> class Descriptor>
int MomentumExchangeBounceBack<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    Dynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue(rho);
    serializer.addValue((int)(fluidDirections.size()));
    serializer.addValues(fluidDirections);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        serializer.addValue(forceIds[iDim]);
    }
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    PLB_PRECONDITION(unserializer.getId() == this->getId());
    Dynamics<T, Descriptor>::unserialize(unserializer);
    unserializer.readValue(rho);
    fluidDirections.resize(unserializer.readValue<int>());
    unserializer.readValues(fluidDirections);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        unserializer.readValue(forceIds[iDim]);
    }
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    enum { numDim = Descriptor<T>::d };
    using namespace indexTemplates;

    // Before collision, the amount of exchanged momentum is computed.
    Array<T, numDim> momentum;
    momentum.resetToZero();
    // Sum over all populations which are incoming from the fluid.
    for (pluint iFluid = 0; iFluid < fluidDirections.size(); ++iFluid) {
        plint iPop = fluidDirections[iFluid];
        plint iOpp = opposite<Descriptor<T> >(iPop);
        for (plint iD = 0; iD < numDim; ++iD) {
            // The momentum contribution is multiplied by two:
            //   One contribution for the momentum loss into the obstacle,
            //   and one contribution for the momentum gain in the subsequent
            //   after-bounce-back streaming step.
            momentum[iD] += (T)2 * Descriptor<T>::c[iPop][iD] * cell[iOpp];
        }
    }
    // Add the momentum exchange for this cell to the total balance due to the
    //   obstacle.
    if (cell.takesStatistics()) {
        for (plint iD = 0; iD < numDim; ++iD) {
            // Add a negative sign, to get the force acting on the obstacle,
            //   and not on the fluid.
            statistics.gatherSum(forceIds[iD], -momentum[iD]);
        }
    }

    // Finally, do the bounce-back operation, which replaces the usual collision.
    for (plint iPop = 1; iPop <= Descriptor<T>::q / 2; ++iPop) {
        std::swap(cell[iPop], cell[iPop + Descriptor<T>::q / 2]);
    }
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    collide(cell, stat);
}

template <typename T, template <typename U> class Descriptor>
T MomentumExchangeBounceBack<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    return T();
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{ }

template <typename T, template <typename U> class Descriptor>
T MomentumExchangeBounceBack<T, Descriptor>::computeDensity(Cell<T, Descriptor> const &cell) const
{
    return rho;
}

template <typename T, template <typename U> class Descriptor>
T MomentumExchangeBounceBack<T, Descriptor>::computePressure(Cell<T, Descriptor> const &cell) const
{
    return T();
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    u.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
T MomentumExchangeBounceBack<T, Descriptor>::computeTemperature(
    Cell<T, Descriptor> const &cell) const
{
    return T();
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::computePiNeq(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    PiNeq.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::computeShearStress(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &stress) const
{
    stress.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::computeHeatFlux(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &q) const
{
    q.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::computeMoment(
    Cell<T, Descriptor> const &cell, plint momentId, T *moment) const
{ }

template <typename T, template <typename U> class Descriptor>
T MomentumExchangeBounceBack<T, Descriptor>::getOmega() const
{
    return T();
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::setOmega(T omega_)
{ }

template <typename T, template <typename U> class Descriptor>
T MomentumExchangeBounceBack<T, Descriptor>::computeRhoBar(Cell<T, Descriptor> const &cell) const
{
    return Descriptor<T>::rhoBar(rho);
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    j.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    rhoBar = Descriptor<T>::rhoBar(rho);
    j.resetToZero();
    PiNeq.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
T MomentumExchangeBounceBack<T, Descriptor>::computeEbar(Cell<T, Descriptor> const &cell) const
{
    return T();
}

template <typename T, template <typename U> class Descriptor>
plint MomentumExchangeBounceBack<T, Descriptor>::numDecomposedVariables(plint order) const
{
    return Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::decompose(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData, plint order) const
{
    rawData.resize(numDecomposedVariables(order));
    std::fill(rawData.begin(), rawData.end(), T());
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::recompose(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData, plint order) const
{ }

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::rescale(
    std::vector<T> &rawData, T xDxInv, T xDt, plint order) const
{ }

template <typename T, template <typename U> class Descriptor>
bool MomentumExchangeBounceBack<T, Descriptor>::hasMoments() const
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
void MomentumExchangeBounceBack<T, Descriptor>::setFluidDirections(
    std::vector<plint> const &fluidDirections_)
{
    fluidDirections = fluidDirections_;
}

template <typename T, template <typename U> class Descriptor>
std::vector<plint> const &MomentumExchangeBounceBack<T, Descriptor>::getFluidDirections() const
{
    return fluidDirections;
}

}  // namespace plb

#endif  // BOUNCE_BACK_MODELS_HH
