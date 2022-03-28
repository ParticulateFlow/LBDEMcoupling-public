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
#ifndef GENERALIZED_BOUNDARY_DYNAMICS_HH
#define GENERALIZED_BOUNDARY_DYNAMICS_HH

#include <Eigen3/Cholesky>
#include <Eigen3/Core>
#include <Eigen3/LU>

#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "generalizedBoundaryDynamics.h"
#include "generalizedBoundaryDynamicsSolvers.h"
#include "generalizedBoundaryDynamicsSolvers.hh"
#include "generalizedCompressibleBoundaryTemplates.h"
#include "generalizedIncompressibleBoundaryTemplates.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

/* *************** Class GeneralizedVelocityBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor>
int GeneralizedVelocityBoundaryDynamics<T, Descriptor>::id = meta::registerGeneralDynamics<
    T, Descriptor, GeneralizedVelocityBoundaryDynamics<T, Descriptor> >(
    "Boundary_GeneralizedVelocity");

template <typename T, template <typename U> class Descriptor>
GeneralizedVelocityBoundaryDynamics<T, Descriptor>::GeneralizedVelocityBoundaryDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
    bool automaticPrepareCollision) :
    StoreVelocityDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision),
    missingIndices(missingIndices_)
{
    knownIndices = indexTemplates::remainingIndexes<Descriptor<T> >(missingIndices);
}

template <typename T, template <typename U> class Descriptor>
GeneralizedVelocityBoundaryDynamics<T, Descriptor>::GeneralizedVelocityBoundaryDynamics(
    HierarchicUnserializer &unserializer) :
    StoreVelocityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GeneralizedVelocityBoundaryDynamics<T, Descriptor>
    *GeneralizedVelocityBoundaryDynamics<T, Descriptor>::clone() const
{
    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GeneralizedVelocityBoundaryDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedVelocityBoundaryDynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    StoreVelocityDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue((int)(missingIndices.size()));
    serializer.addValues(missingIndices);
    serializer.addValue((int)(knownIndices.size()));
    serializer.addValues(knownIndices);
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedVelocityBoundaryDynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    PLB_PRECONDITION(unserializer.getId() == this->getId());
    StoreVelocityDynamics<T, Descriptor>::unserialize(unserializer);
    missingIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(missingIndices);
    knownIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(knownIndices);
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedVelocityBoundaryDynamics<T, Descriptor>::computeUlb(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &uLb) const
{
    this->computeVelocity(cell, uLb);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        uLb[iDim] -= 0.5 * getExternalForceComponent(cell, iDim);
    }
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedVelocityBoundaryDynamics<T, Descriptor>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    Array<T, Descriptor<T>::d> uLb;
    computeUlb(cell, uLb);

    DirichletVelocityBoundarySolver<T, Descriptor> bc(missingIndices, knownIndices, uLb);
    bc.apply(cell, this->getBaseDynamics());
}

/* *************** Class GeneralizedMassConservingVelocityBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor>
int GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor> >(
        "Boundary_MassConservingGeneralizedVelocity");

template <typename T, template <typename U> class Descriptor>
GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>::
    GeneralizedMassConservingVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        std::vector<plint> knownIndices_, std::vector<plint> inGoingIndices_,
        bool automaticPrepareCollision) :
    StoreVelocityDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision),
    missingIndices(missingIndices_),
    knownIndices(knownIndices_),
    inGoingIndices(inGoingIndices_)
{ }

template <typename T, template <typename U> class Descriptor>
GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>::
    GeneralizedMassConservingVelocityBoundaryDynamics(HierarchicUnserializer &unserializer) :
    StoreVelocityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>
    *GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>::clone() const
{
    return new GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    StoreVelocityDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue((int)(missingIndices.size()));
    serializer.addValues(missingIndices);
    serializer.addValue((int)(knownIndices.size()));
    serializer.addValues(knownIndices);
    serializer.addValue((int)(inGoingIndices.size()));
    serializer.addValues(inGoingIndices);
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    PLB_PRECONDITION(unserializer.getId() == this->getId());
    StoreVelocityDynamics<T, Descriptor>::unserialize(unserializer);
    missingIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(missingIndices);
    knownIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(knownIndices);
    inGoingIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(inGoingIndices);
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>::computeUlb(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &uLb) const
{
    this->computeVelocity(cell, uLb);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        uLb[iDim] -= 0.5 * getExternalForceComponent(cell, iDim);
    }
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    Array<T, Descriptor<T>::d> uLb;
    computeUlb(cell, uLb);

    DirichletMassConservingVelocityBoundarySolver<T, Descriptor> bc(
        missingIndices, knownIndices, inGoingIndices, uLb);
    bc.apply(cell, this->getBaseDynamics());
}

/* *************** Class GeneralizedDensityBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_GeneralizedDensity_") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::
    GeneralizedDensityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        bool automaticPrepareCollision_) :
    DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics_, automaticPrepareCollision_),
    missingIndices(missingIndices_)
{
    knownIndices = indexTemplates::remainingIndexes<Descriptor<T> >(missingIndices);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::
    GeneralizedDensityBoundaryDynamics(HierarchicUnserializer &unserializer) :
    DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>
    *GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::serialize(
    HierarchicSerializer &serializer) const
{
    DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::serialize(serializer);
    serializer.addValue((int)(missingIndices.size()));
    serializer.addValues(missingIndices);
    serializer.addValue((int)(knownIndices.size()));
    serializer.addValues(knownIndices);
    //     serializer.addValues(u);
    //     serializer.addValues(PiNeq);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::unserialize(
    HierarchicUnserializer &unserializer)
{
    PLB_PRECONDITION(unserializer.getId() == this->getId());
    DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation>::unserialize(
        unserializer);
    missingIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(missingIndices);
    knownIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(knownIndices);
    //     unserializer.readValues(u);
    //     unserializer.readValues(PiNeq);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> u, j;
    this->computeRhoBarJ(cell, rhoBar, j);
    T rho = Descriptor<T>::fullRho(rhoBar);
    T invRho = Descriptor<T>::invRho(rhoBar);
    u = invRho * j;

    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    PiNeq.resetToZero();

    DirichletDensityBoundarySolver<T, Descriptor, direction> bc(
        missingIndices, knownIndices, rho, u, PiNeq, 1.0e-7);
    bc.apply(cell, this->getBaseDynamics());
}

/* ********** Class GeneralizedVelocityTemperatureBoundaryDynamics ********* */

template <typename T, template <typename U> class Descriptor>
int GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor> >(
        "Boundary_GeneralizedVelocityTemperature");

template <typename T, template <typename U> class Descriptor>
GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>::
    GeneralizedVelocityTemperatureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        bool massConserving_, bool automaticPrepareCollision) :
    StoreTemperatureAndVelocityDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision),
    missingIndices(missingIndices_),
    massConserving(massConserving_)
{
    knownIndices = indexTemplates::remainingIndexes<Descriptor<T> >(missingIndices);
}

template <typename T, template <typename U> class Descriptor>
GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>::
    GeneralizedVelocityTemperatureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
        std::vector<plint> knownIndices_, bool massConserving_, bool automaticPrepareCollision) :
    StoreTemperatureAndVelocityDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision),
    missingIndices(missingIndices_),
    knownIndices(knownIndices_),
    massConserving(massConserving_)
{ }

template <typename T, template <typename U> class Descriptor>
GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>::
    GeneralizedVelocityTemperatureBoundaryDynamics(HierarchicUnserializer &unserializer) :
    StoreTemperatureAndVelocityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>
    *GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>::clone() const
{
    return new GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    StoreTemperatureAndVelocityDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue((int)(missingIndices.size()));
    serializer.addValues(missingIndices);
    serializer.addValue((int)(knownIndices.size()));
    serializer.addValues(knownIndices);
    serializer.addValue(massConserving);
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    PLB_PRECONDITION(unserializer.getId() == this->getId());
    StoreTemperatureAndVelocityDynamics<T, Descriptor>::unserialize(unserializer);
    missingIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(missingIndices);
    knownIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(knownIndices);
    massConserving = unserializer.readValue<bool>();
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>::computeUlb(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &uLb) const
{
    this->computeVelocity(cell, uLb);
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        uLb[iDim] -= 0.5 * getExternalForceComponent(cell, iDim);
    }
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedVelocityTemperatureBoundaryDynamics<T, Descriptor>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    Array<T, Descriptor<T>::d> uLb;
    computeUlb(cell, uLb);
    T thetaBar = this->computeTemperature(cell) * Descriptor<T>::invCs2 - (T)1;

    T rho;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    Array<T, SymmetricRankThreeTensor<T, Descriptor>::n> qNeq;
    //     generalizedComprTempBoundaryTemplates<T,Descriptor>::
    //         solveLinearSystemTrLessPiNeq(cell, uLb, thetaBar,missingIndices, knownIndices, rho,
    //         PiNeq, qNeq, massConserving);

    generalizedComprTempBoundaryTemplates<T, Descriptor>::solveLinearSystem(
        cell, uLb, thetaBar, missingIndices, knownIndices, rho, PiNeq, qNeq, massConserving);

    T rhoBar = Descriptor<T>::rhoBar(rho);

    Array<T, Descriptor<T>::d> j;
    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        j[iDim] = rho * uLb[iDim];
    }
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    // TODO not elegant. Regularize functions should be made more generic
    typedef Descriptor<T> L;
    Array<T, Descriptor<T>::q> f, feq;
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        cell[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
        cell[iPop] += offEquilibriumTemplates<T, Descriptor>::fromPiAndQtoFneq(iPop, PiNeq, qNeq);
    }
}

/* *************** Class GeneralizedNextToWallBoundaryDynamics ************* */

template <typename T, template <typename U> class Descriptor>
int GeneralizedNextToBoundaryDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, GeneralizedNextToBoundaryDynamics<T, Descriptor> >(
        "Boundary_GeneralizedNextToBoundaryDynamics");

template <typename T, template <typename U> class Descriptor>
GeneralizedNextToBoundaryDynamics<T, Descriptor>::GeneralizedNextToBoundaryDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
    bool automaticPrepareCollision) :
    BoundaryCompositeDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision),
    missingIndices(missingIndices_)
{
    knownIndices = indexTemplates::remainingIndexes<Descriptor<T> >(missingIndices);

    rho = (T)1;
    thetaBar = T();
    u.resetToZero();
    PiNeq.resetToZero();
    Qneq.resetToZero();
}

template <typename T, template <typename U> class Descriptor>
GeneralizedNextToBoundaryDynamics<T, Descriptor>::GeneralizedNextToBoundaryDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, std::vector<plint> missingIndices_,
    std::vector<plint> knownIndices_, bool automaticPrepareCollision) :
    BoundaryCompositeDynamics<T, Descriptor>(baseDynamics_, automaticPrepareCollision),
    missingIndices(missingIndices_),
    knownIndices(knownIndices_)
{ }

template <typename T, template <typename U> class Descriptor>
GeneralizedNextToBoundaryDynamics<T, Descriptor>::GeneralizedNextToBoundaryDynamics(
    HierarchicUnserializer &unserializer) :
    BoundaryCompositeDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
GeneralizedNextToBoundaryDynamics<T, Descriptor>
    *GeneralizedNextToBoundaryDynamics<T, Descriptor>::clone() const
{
    return new GeneralizedNextToBoundaryDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int GeneralizedNextToBoundaryDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedNextToBoundaryDynamics<T, Descriptor>::serialize(
    HierarchicSerializer &serializer) const
{
    BoundaryCompositeDynamics<T, Descriptor>::serialize(serializer);
    serializer.addValue((int)(missingIndices.size()));
    serializer.addValues(missingIndices);
    serializer.addValue((int)(knownIndices.size()));
    serializer.addValues(knownIndices);
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedNextToBoundaryDynamics<T, Descriptor>::unserialize(
    HierarchicUnserializer &unserializer)
{
    PLB_PRECONDITION(unserializer.getId() == this->getId());
    BoundaryCompositeDynamics<T, Descriptor>::unserialize(unserializer);
    missingIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(missingIndices);
    knownIndices.resize(unserializer.readValue<int>());
    unserializer.readValues(knownIndices);
}

template <typename T, template <typename U> class Descriptor>
void GeneralizedNextToBoundaryDynamics<T, Descriptor>::completePopulations(
    Cell<T, Descriptor> &cell)
{
    T epsilon = 1.0e-10;
    T resSum;
    bool converged = generalizedComprTempBoundaryTemplates<T, Descriptor>::iterativelySolveSystem(
        cell, rho, u, thetaBar, PiNeq, Qneq, knownIndices, missingIndices, epsilon, resSum);
    //     bool converged = generalizedComprTempBoundaryTemplates<T,Descriptor>::
    //         iterativelySolveSystemTrLessPiNeq(cell, rho, u, thetaBar, PiNeq, Qneq, knownIndices,
    //         missingIndices, epsilon, resSum);
    if (!converged) {
        pcout << "Never converged... Exiting." << std::endl;
        exit(-1);
    }
    Array<T, Descriptor<T>::d> j;
    j = rho * u;
    T rhoBar = Descriptor<T>::rhoBar(rho);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
        cell[iPop] += offEquilibriumTemplates<T, Descriptor>::fromPiAndQtoFneq(iPop, PiNeq, Qneq);
    }
}

}  // namespace plb

#endif  // GENERALIZED_BOUNDARY_DYNAMICS_HH
