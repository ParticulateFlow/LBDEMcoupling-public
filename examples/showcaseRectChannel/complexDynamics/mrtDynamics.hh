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

/* Orestis Malaspinas contributed this code. */

/** \file
 * MRT dynamics -- generic implementation.
 */
#ifndef MRT_DYNAMICS_HH
#define MRT_DYNAMICS_HH

#include <algorithm>
#include <limits>

#include "complexDynamics/mrtDynamics.h"
#include "core/latticeStatistics.h"
#include "latticeBoltzmann/d3q13Templates.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/mrtTemplates.h"

namespace plb {

/* *************** Class MRTdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int MRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, MRTdynamics<T, Descriptor> >("MRT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
MRTdynamics<T, Descriptor>::MRTdynamics(T omega_) : IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
MRTdynamics<T, Descriptor>::MRTdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
MRTdynamics<T, Descriptor> *MRTdynamics<T, Descriptor>::clone() const
{
    return new MRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int MRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void MRTdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T jSqr = mrtTemp::mrtCollision(cell, this->getOmega());

    if (cell.takesStatistics()) {
        T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
        gatherStatistics(
            statistics, rhoBar,
            jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar));
    }
}

template <typename T, template <typename U> class Descriptor>
void MRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T jSqr = mrtTemp::mrtCollision(cell, rhoBar, j, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(
            statistics, rhoBar,
            jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar));
    }
}

template <typename T, template <typename U> class Descriptor>
T MRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class IncMRTdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int IncMRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, IncMRTdynamics<T, Descriptor> >("IncMRT");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
IncMRTdynamics<T, Descriptor>::IncMRTdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
IncMRTdynamics<T, Descriptor>::IncMRTdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
IncMRTdynamics<T, Descriptor> *IncMRTdynamics<T, Descriptor>::clone() const
{
    return new IncMRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int IncMRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void IncMRTdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T jSqr = mrtTemp::incMrtCollision(cell, this->getOmega());

    if (cell.takesStatistics()) {
        T rhoBar = momentTemplates<T, Descriptor>::get_rhoBar(cell);
        gatherStatistics(
            statistics, rhoBar,
            jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar));
    }
}

template <typename T, template <typename U> class Descriptor>
void IncMRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    typedef mrtTemplates<T, Descriptor> mrtTemp;

    T jSqr = mrtTemp::incMrtCollision(cell, rhoBar, j, this->getOmega());

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T IncMRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
bool IncMRTdynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

template <typename T, template <typename U> class Descriptor>
void IncMRTdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T dummyRhoBar;
    this->computeRhoBarJ(cell, dummyRhoBar, u);
}

template <typename T, template <typename U> class Descriptor>
void IncMRTdynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    // Incompressible: rho0=1
    T invRho0 = (T)1.;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq, invRho0);
}

/* *************** Class D3Q13Dynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int D3Q13Dynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, D3Q13Dynamics<T, Descriptor> >("D3Q13");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
D3Q13Dynamics<T, Descriptor>::D3Q13Dynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
D3Q13Dynamics<T, Descriptor>::D3Q13Dynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
D3Q13Dynamics<T, Descriptor> *D3Q13Dynamics<T, Descriptor>::clone() const
{
    return new D3Q13Dynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int D3Q13Dynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void D3Q13Dynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef d3q13Templates<T> d3q13Temp;
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T rho = Descriptor<T>::fullRho(rhoBar);
    Array<T, Descriptor<T>::d> u(j / rho);

    T omega = this->getOmega();
    T lambda_nu = (T)2 / ((T)8 / omega - (T)3);
    T lambda_nu_prime = (T)2 / ((T)4 / omega - (T)1);

    T uSqr = d3q13Temp::collision(cell, rho, u, lambda_nu, lambda_nu_prime);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void D3Q13Dynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    typedef d3q13Templates<T> d3q13Temp;
    T rho = Descriptor<T>::fullRho(rhoBar);
    Array<T, Descriptor<T>::d> u(j / rho);

    T omega = this->getOmega();
    T lambda_nu = (T)2 / ((T)8 / omega - (T)3);
    T lambda_nu_prime = (T)2 / ((T)4 / omega - (T)1);

    T uSqr = d3q13Temp::collision(cell, rho, u, lambda_nu, lambda_nu_prime);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T D3Q13Dynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
bool D3Q13Dynamics<T, Descriptor>::velIsJ() const
{
    return false;
}

}  // namespace plb

#endif  // MRT_DYNAMICS_HH
