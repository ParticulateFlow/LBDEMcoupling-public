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
 * Parallel dynamics object -- generic template code
 */
#ifndef PARALLEL_DYNAMICS_HH
#define PARALLEL_DYNAMICS_HH

#include "parallelism/mpiManager.h"
#include "parallelism/parallelDynamics.h"

namespace plb {

#ifdef PLB_MPI_PARALLEL

////////////////////// Class ParallelDynamics /////////////////////////////

template <typename T, template <typename U> class Descriptor>
ParallelDynamics<T, Descriptor>::ParallelDynamics(
    std::vector<Cell<T, Descriptor> *> &baseCells_, bool hasBulkCell_) :
    baseCells(baseCells_), hasBulkCell(hasBulkCell_)
{ }

template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> *ParallelDynamics<T, Descriptor>::clone() const
{
    return new ParallelDynamics(*this);
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics_)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->collide(statistics_);
    }
}

template <typename T, template <typename U> class Descriptor>
T ParallelDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T eq = T();
    if (hasBulkCell) {
        eq = baseCells[0]->computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
    }
    global::mpi().bCastThroughMaster(&eq, 1, hasBulkCell);
    return eq;
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->regularize(rhoBar, j, jSqr, PiNeq, thetaBar);
    }
}

template <typename T, template <typename U> class Descriptor>
T ParallelDynamics<T, Descriptor>::computeDensity(Cell<T, Descriptor> const &cell) const
{
    T rho = T();
    if (hasBulkCell) {
        rho = baseCells[0]->computeDensity();
    }
    global::mpi().bCastThroughMaster(&rho, 1, hasBulkCell);
    return rho;
}

template <typename T, template <typename U> class Descriptor>
T ParallelDynamics<T, Descriptor>::computePressure(Cell<T, Descriptor> const &cell) const
{
    T p = T();
    if (hasBulkCell) {
        p = baseCells[0]->computePressure();
    }
    global::mpi().bCastThroughMaster(&p, 1, hasBulkCell);
    return p;
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    if (hasBulkCell) {
        baseCells[0]->computeVelocity(u);
    }
    global::mpi().bCastThroughMaster(&u[0], Descriptor<T>::d, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
T ParallelDynamics<T, Descriptor>::computeTemperature(Cell<T, Descriptor> const &cell) const
{
    T theta = T();
    if (hasBulkCell) {
        theta = baseCells[0]->computeTemperature();
    }
    global::mpi().bCastThroughMaster(&theta, 1, hasBulkCell);
    return theta;
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::computePiNeq(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().computePiNeq(*baseCells[0], PiNeq);
    }
    global::mpi().bCastThroughMaster(&PiNeq[0], SymmetricTensor<T, Descriptor>::n, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::computeShearStress(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &stress) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().computeShearStress(*baseCells[0], stress);
    }
    global::mpi().bCastThroughMaster(&stress[0], SymmetricTensor<T, Descriptor>::n, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::computeHeatFlux(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &q) const
{
    if (hasBulkCell) {
        baseCells[0]->computeHeatFlux(q);
    }
    global::mpi().bCastThroughMaster(&q[0], Descriptor<T>::d, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::computeMoment(
    Cell<T, Descriptor> const &cell, plint momentId, T *moment) const
{
    // Cannot transfer generic moment through MPI,
    // because the type size is unknown
    PLB_PRECONDITION(false);
}

template <typename T, template <typename U> class Descriptor>
T ParallelDynamics<T, Descriptor>::getOmega() const
{
    T omega = T();
    if (hasBulkCell) {
        omega = baseCells[0]->getDynamics().getOmega();
    }
    global::mpi().bCastThroughMaster(&omega, 1, hasBulkCell);
    return omega;
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::setOmega(T omega_)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->getDynamics().setOmega(omega_);
    }
}

template <typename T, template <typename U> class Descriptor>
T ParallelDynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    T parameter = T();
    if (hasBulkCell) {
        parameter = baseCells[0]->getDynamics().getParameter(whichParameter);
    }
    global::mpi().bCastThroughMaster(&parameter, 1, hasBulkCell);
    return parameter;
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::setParameter(plint whichParameter, T value)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->getDynamics().setParameter(whichParameter, value);
    }
}

template <typename T, template <typename U> class Descriptor>
plint ParallelDynamics<T, Descriptor>::numDecomposedVariables(plint order) const
{
    plint numVariables = 0;
    if (hasBulkCell) {
        numVariables = baseCells[0]->getDynamics().numDecomposedVariables(order);
    }
    global::mpi().bCastThroughMaster(&numVariables, 1, hasBulkCell);
    return numVariables;
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::decompose(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData, plint order) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().decompose(cell, rawData, order);
    }
    plint numVariables = numDecomposedVariables(order);
    global::mpi().bCastThroughMaster(&rawData[0], numVariables, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::recompose(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData, plint order) const
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->getDynamics().recompose(cell, rawData, order);
    }
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::rescale(
    std::vector<T> &rawData, T xDxInv, T xDt, plint order) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().rescale(rawData, xDxInv, xDt, order);
    }
    plint numVariables = numDecomposedVariables(order);
    global::mpi().bCastThroughMaster(&rawData[0], numVariables, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::getPopulations(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::q> &f) const
{
    if (hasBulkCell) {
        baseCells[0]->getPopulations(f);
    }
    global::mpi().bCastThroughMaster(&f[0], Descriptor<T>::q, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::getExternalField(
    Cell<T, Descriptor> const &cell, plint pos, plint size, T *ext) const
{
    if (hasBulkCell) {
        baseCells[0]->getExternalField(pos, size, ext);
    }
    global::mpi().bCastThroughMaster(ext, Descriptor<T>::ExternalField::numScalars, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::setPopulations(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::q> const &f)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->setPopulations(f);
    }
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::setExternalField(
    Cell<T, Descriptor> &cell, plint pos, plint size, const T *ext)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->setExternalField(pos, size, ext);
    }
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::defineDensity(Cell<T, Descriptor> &cell, T rho)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->defineDensity(rho);
    }
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::defineVelocity(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &u)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->defineVelocity(u);
    }
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::defineTemperature(Cell<T, Descriptor> &cell, T temperature)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->defineTemperature(temperature);
    }
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::defineHeatFlux(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &q)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->defineHeatFlux(q);
    }
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::definePiNeq(
    Cell<T, Descriptor> &cell, Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->definePiNeq(PiNeq);
    }
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::defineMoment(
    Cell<T, Descriptor> &cell, plint momentId, T const *value)
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->defineMoment(momentId, value);
    }
}

template <typename T, template <typename U> class Descriptor>
T ParallelDynamics<T, Descriptor>::computeRhoBar(Cell<T, Descriptor> const &cell) const
{
    T rhoBar = T();
    if (hasBulkCell) {
        baseCells[0]->getDynamics().computeRhoBar(*baseCells[0]);
    }
    global::mpi().bCastThroughMaster(&rhoBar, 1, hasBulkCell);
    return rhoBar;
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().computeRhoBarJ(*baseCells[0], rhoBar, j);
    }
    global::mpi().bCastThroughMaster(&rhoBar, 1, hasBulkCell);
    global::mpi().bCastThroughMaster(&j[0], Descriptor<T>::d, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ParallelDynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().computeRhoBarJPiNeq(*baseCells[0], rhoBar, j, PiNeq);
    }
    global::mpi().bCastThroughMaster(&rhoBar, 1, hasBulkCell);
    global::mpi().bCastThroughMaster(&j[0], Descriptor<T>::d, hasBulkCell);
    global::mpi().bCastThroughMaster(&PiNeq[0], SymmetricTensor<T, Descriptor>::n, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
T ParallelDynamics<T, Descriptor>::computeEbar(Cell<T, Descriptor> const &cell) const
{
    T eBar = T();
    if (hasBulkCell) {
        eBar = baseCells[0]->getDynamics().computeEbar(*baseCells[0]);
    }
    global::mpi().bCastThroughMaster(&eBar, 1, hasBulkCell);
    return eBar;
}

////////////////// Class ConstParallelDynamics /////////////////////////

template <typename T, template <typename U> class Descriptor>
ConstParallelDynamics<T, Descriptor>::ConstParallelDynamics(
    std::vector<Cell<T, Descriptor> const *> &baseCells_, bool hasBulkCell_) :
    baseCells(baseCells_), hasBulkCell(hasBulkCell_)
{ }

template <typename T, template <typename U> class Descriptor>
Dynamics<T, Descriptor> *ConstParallelDynamics<T, Descriptor>::clone() const
{
    return new ConstParallelDynamics(*this);
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics_)
{ }

template <typename T, template <typename U> class Descriptor>
T ConstParallelDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T eq = T();
    if (hasBulkCell) {
        eq = baseCells[0]->computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
    }
    global::mpi().bCastThroughMaster(&eq, 1, hasBulkCell);
    return eq;
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::regularize(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr,
    Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq, T thetaBar) const
{ }

template <typename T, template <typename U> class Descriptor>
T ConstParallelDynamics<T, Descriptor>::computeDensity(Cell<T, Descriptor> const &cell) const
{
    T rho = T();
    if (hasBulkCell) {
        rho = baseCells[0]->computeDensity();
    }
    global::mpi().bCastThroughMaster(&rho, 1, hasBulkCell);
    return rho;
}

template <typename T, template <typename U> class Descriptor>
T ConstParallelDynamics<T, Descriptor>::computePressure(Cell<T, Descriptor> const &cell) const
{
    T p = T();
    if (hasBulkCell) {
        p = baseCells[0]->computePressure();
    }
    global::mpi().bCastThroughMaster(&p, 1, hasBulkCell);
    return p;
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    if (hasBulkCell) {
        baseCells[0]->computeVelocity(u);
    }
    global::mpi().bCastThroughMaster(&u[0], Descriptor<T>::d, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
T ConstParallelDynamics<T, Descriptor>::computeTemperature(Cell<T, Descriptor> const &cell) const
{
    T theta = T();
    if (hasBulkCell) {
        theta = baseCells[0]->computeTemperature();
    }
    global::mpi().bCastThroughMaster(&theta, 1, hasBulkCell);
    return theta;
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::computePiNeq(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().computePiNeq(*baseCells[0], PiNeq);
    }
    global::mpi().bCastThroughMaster(&PiNeq[0], SymmetricTensor<T, Descriptor>::n, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::computeShearStress(
    Cell<T, Descriptor> const &cell, Array<T, SymmetricTensor<T, Descriptor>::n> &stress) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().computeShearStress(*baseCells[0], stress);
    }
    global::mpi().bCastThroughMaster(&stress[0], SymmetricTensor<T, Descriptor>::n, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::computeHeatFlux(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &q) const
{
    if (hasBulkCell) {
        baseCells[0]->computeHeatFlux(q);
    }
    global::mpi().bCastThroughMaster(&q[0], Descriptor<T>::d, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::computeMoment(
    Cell<T, Descriptor> const &cell, plint momentId, T *moment) const
{
    // Cannot transfer generic moment through MPI,
    // because the type size is unknown
    PLB_PRECONDITION(false);
}

template <typename T, template <typename U> class Descriptor>
T ConstParallelDynamics<T, Descriptor>::getOmega() const
{
    T omega = T();
    if (hasBulkCell) {
        omega = baseCells[0]->getDynamics().getOmega();
    }
    global::mpi().bCastThroughMaster(&omega, 1, hasBulkCell);
    return omega;
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::setOmega(T omega_)
{ }

template <typename T, template <typename U> class Descriptor>
T ConstParallelDynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    T parameter = T();
    if (hasBulkCell) {
        parameter = baseCells[0]->getDynamics().getParameter(whichParameter);
    }
    global::mpi().bCastThroughMaster(&parameter, 1, hasBulkCell);
    return parameter;
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::setParameter(plint whichParameter, T value)
{ }

template <typename T, template <typename U> class Descriptor>
plint ConstParallelDynamics<T, Descriptor>::numDecomposedVariables(plint order) const
{
    plint numVariables = 0;
    if (hasBulkCell) {
        numVariables = baseCells[0]->getDynamics().numDecomposedVariables(order);
    }
    global::mpi().bCastThroughMaster(&numVariables, 1, hasBulkCell);
    return numVariables;
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::decompose(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData, plint order) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().decompose(cell, rawData, order);
    }
    plint numVariables = numDecomposedVariables(order);
    global::mpi().bCastThroughMaster(&rawData[0], numVariables, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::recompose(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData, plint order) const
{
    for (pluint iCell = 0; iCell < baseCells.size(); ++iCell) {
        baseCells[iCell]->getDynamics().recompose(cell, rawData, order);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::rescale(
    std::vector<T> &rawData, T xDxInv, T xDt, plint order) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().rescale(rawData, xDxInv, xDt, order);
    }
    plint numVariables = numDecomposedVariables(order);
    global::mpi().bCastThroughMaster(&rawData[0], numVariables, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::getPopulations(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::q> &f) const
{
    if (hasBulkCell) {
        baseCells[0]->getPopulations(f);
    }
    global::mpi().bCastThroughMaster(&f[0], Descriptor<T>::q, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::getExternalField(
    Cell<T, Descriptor> const &cell, plint pos, plint size, T *ext) const
{
    if (hasBulkCell) {
        baseCells[0]->getExternalField(pos, size, ext);
    }
    global::mpi().bCastThroughMaster(ext, Descriptor<T>::ExternalField::numScalars, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::setPopulations(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::q> const &f)
{ }

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::setExternalField(
    Cell<T, Descriptor> &cell, plint pos, plint size, const T *ext)
{ }

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::defineDensity(Cell<T, Descriptor> &cell, T rho)
{ }

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::defineVelocity(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &u)
{ }

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::defineTemperature(
    Cell<T, Descriptor> &cell, T temperature)
{ }

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::defineHeatFlux(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::d> const &q)
{ }

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::definePiNeq(
    Cell<T, Descriptor> &cell, Array<T, SymmetricTensor<T, Descriptor>::n> const &PiNeq)
{ }

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::defineMoment(
    Cell<T, Descriptor> &cell, plint momentId, T const *value)
{ }

template <typename T, template <typename U> class Descriptor>
T ConstParallelDynamics<T, Descriptor>::computeRhoBar(Cell<T, Descriptor> const &cell) const
{
    T rhoBar = T();
    if (hasBulkCell) {
        baseCells[0]->getDynamics().computeRhoBar(*baseCells[0]);
    }
    global::mpi().bCastThroughMaster(&rhoBar, 1, hasBulkCell);
    return rhoBar;
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::computeRhoBarJ(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().computeRhoBarJ(*baseCells[0], rhoBar, j);
    }
    global::mpi().bCastThroughMaster(&rhoBar, 1, hasBulkCell);
    global::mpi().bCastThroughMaster(&j[0], Descriptor<T>::d, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
void ConstParallelDynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    if (hasBulkCell) {
        baseCells[0]->getDynamics().computeRhoBarJPiNeq(*baseCells[0], rhoBar, j, PiNeq);
    }
    global::mpi().bCastThroughMaster(&rhoBar, 1, hasBulkCell);
    global::mpi().bCastThroughMaster(&j[0], Descriptor<T>::d, hasBulkCell);
    global::mpi().bCastThroughMaster(&PiNeq[0], SymmetricTensor<T, Descriptor>::n, hasBulkCell);
}

template <typename T, template <typename U> class Descriptor>
T ConstParallelDynamics<T, Descriptor>::computeEbar(Cell<T, Descriptor> const &cell) const
{
    T eBar = T();
    if (hasBulkCell) {
        eBar = baseCells[0]->getDynamics().computeEbar(*baseCells[0]);
    }
    global::mpi().bCastThroughMaster(&eBar, 1, hasBulkCell);
    return eBar;
}

#endif

}  // namespace plb

#endif  // defined MULTIBLOCK_DYNAMICS_H
