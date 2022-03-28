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
 * Implementation of the LW-ACM model by Pietro Asinari and others -- generic implementation.
 */
#ifndef ASINARI_MODEL_HH
#define ASINARI_MODEL_HH

#include <algorithm>
#include <limits>

#include "complexDynamics/asinariModel.h"
#include "complexDynamics/asinariTemplates.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"
#include "finiteDifference/fdStencils1D.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"

namespace plb {

/* *************** Class AsinariDynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int AsinariDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, AsinariDynamics<T, Descriptor> >("Asinari");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
AsinariDynamics<T, Descriptor>::AsinariDynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{
    computePrefactor();
}

template <typename T, template <typename U> class Descriptor>
AsinariDynamics<T, Descriptor>::AsinariDynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
AsinariDynamics<T, Descriptor> *AsinariDynamics<T, Descriptor>::clone() const
{
    return new AsinariDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int AsinariDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void AsinariDynamics<T, Descriptor>::setOmega(T omega_)
{
    IsoThermalBulkDynamics<T, Descriptor>::setOmega(omega_);
    computePrefactor();
}

template <typename T, template <typename U> class Descriptor>
void AsinariDynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef Descriptor<T> D;
    static const int rhoBarOfs = D::ExternalField::rhoBarBeginsAt;
    static const int jOfs = D::ExternalField::jBeginsAt;
    typedef asinariTemplates<T, Descriptor> asi;
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);

    *(cell.getExternal(rhoBarOfs)) = rhoBar;
    j.to_cArray(cell.getExternal(jOfs));
    T invRho = D::invRho(rhoBar);
    T jSqr = normSqr(j);

    T uSqr = asi::bgk_collision_stage1(cell, rhoBar, invRho, j, jSqr, prefactor);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void AsinariDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    typedef Descriptor<T> D;
    static const int rhoBarOfs = D::ExternalField::rhoBarBeginsAt;
    static const int jOfs = D::ExternalField::jBeginsAt;
    typedef asinariTemplates<T, Descriptor> asi;

    *(cell.getExternal(rhoBarOfs)) = rhoBar;
    j.to_cArray(cell.getExternal(jOfs));
    T invRho = D::invRho(rhoBar);
    T jSqr = normSqr(j);

    T uSqr = asi::bgk_collision_stage1(cell, rhoBar, invRho, j, jSqr, prefactor);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T AsinariDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void AsinariDynamics<T, Descriptor>::computePrefactor()
{
    T omega = this->getOmega();
    prefactor = (T)2 * (omega - (T)1) / omega;
}

template <typename T, template <typename U> class Descriptor>
T AsinariDynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    if (whichParameter == 1000) {
        return prefactor;
    } else {
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

/* *************** Class IncAsinariDynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int IncAsinariDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, IncAsinariDynamics<T, Descriptor> >(
        "Asinari_Incompressible");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
IncAsinariDynamics<T, Descriptor>::IncAsinariDynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{
    computePrefactor();
}

template <typename T, template <typename U> class Descriptor>
IncAsinariDynamics<T, Descriptor>::IncAsinariDynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
IncAsinariDynamics<T, Descriptor> *IncAsinariDynamics<T, Descriptor>::clone() const
{
    return new IncAsinariDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int IncAsinariDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void IncAsinariDynamics<T, Descriptor>::setOmega(T omega_)
{
    IsoThermalBulkDynamics<T, Descriptor>::setOmega(omega_);
    computePrefactor();
}

template <typename T, template <typename U> class Descriptor>
void IncAsinariDynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T dummyRhoBar;
    this->computeRhoBarJ(cell, dummyRhoBar, u);
}

template <typename T, template <typename U> class Descriptor>
void IncAsinariDynamics<T, Descriptor>::computeRhoBarJPiNeq(
    Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j,
    Array<T, SymmetricTensor<T, Descriptor>::n> &PiNeq) const
{
    // Incompressible: rho0=1.
    T invRho0 = (T)1.0;
    momentTemplates<T, Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq, invRho0);
}

template <typename T, template <typename U> class Descriptor>
void IncAsinariDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef Descriptor<T> D;
    static const int rhoBarOfs = D::ExternalField::rhoBarBeginsAt;
    static const int jOfs = D::ExternalField::jBeginsAt;
    typedef asinariTemplates<T, Descriptor> asi;
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);

    *(cell.getExternal(rhoBarOfs)) = rhoBar;
    j.to_cArray(cell.getExternal(jOfs));
    // For the incompressible BGK dynamics, the "1/rho" pre-factor of
    // the O(Ma^2) term is unity.
    T invRho = (T)1;
    T jSqr = normSqr(j);

    T uSqr = asi::bgk_collision_stage1(cell, rhoBar, invRho, j, jSqr, prefactor);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void IncAsinariDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    typedef Descriptor<T> D;
    static const int rhoBarOfs = D::ExternalField::rhoBarBeginsAt;
    static const int jOfs = D::ExternalField::jBeginsAt;
    typedef asinariTemplates<T, Descriptor> asi;

    *(cell.getExternal(rhoBarOfs)) = rhoBar;
    j.to_cArray(cell.getExternal(jOfs));
    // For the incompressible BGK dynamics, the "1/rho" pre-factor of
    // the O(Ma^2) term is unity.
    T invRho = (T)1;
    T jSqr = normSqr(j);

    T uSqr = asi::bgk_collision_stage1(cell, rhoBar, invRho, j, jSqr, prefactor);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T IncAsinariDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    // For the incompressible BGK dynamics, the "1/rho" pre-factor of
    // the O(Ma^2) term is unity.
    T invRho = (T)1;
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void IncAsinariDynamics<T, Descriptor>::computePrefactor()
{
    T omega = this->getOmega();
    prefactor = (T)2 * (omega - (T)1) / omega;
}

template <typename T, template <typename U> class Descriptor>
T IncAsinariDynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    if (whichParameter == 1000) {
        return prefactor;
    } else {
        return IsoThermalBulkDynamics<T, Descriptor>::getParameter(whichParameter);
    }
}

template <typename T, template <typename U> class Descriptor>
bool IncAsinariDynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

/* ************* Class AsinariPostCollide3D ******************* */

template <typename T, template <typename U> class Descriptor>
void AsinariPostCollide3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    typedef Descriptor<T> D;
    typedef asinariTemplates<T, Descriptor> asi;
    static const int jOfs = D::ExternalField::jBeginsAt;
    static const int bounceBackID = BounceBack<T, Descriptor>().getId();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                if (cell.getDynamics().getId() != bounceBackID) {
                    T prefactor = cell.getDynamics().getParameter(1000);
                    Array<T, 3> j;
                    j.from_cArray(cell.getExternal(jOfs));
                    asi::bgk_collision_stage3(cell, j, prefactor);
                }
            }
        }
    }
}

/* ************* Class AsinariPostCollide2D ******************* */

template <typename T, template <typename U> class Descriptor>
void AsinariPostCollide2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    typedef Descriptor<T> D;
    typedef asinariTemplates<T, Descriptor> asi;
    static const int jOfs = D::ExternalField::jBeginsAt;
    static const int bounceBackID = BounceBack<T, Descriptor>().getId();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> &cell = lattice.get(iX, iY);
            if (cell.getDynamics().getId() != bounceBackID) {
                T prefactor = cell.getDynamics().getParameter(1000);
                Array<T, 2> j;
                j.from_cArray(cell.getExternal(jOfs));
                asi::bgk_collision_stage3(cell, j, prefactor);
            }
        }
    }
}

}  // namespace plb

#endif  // ASINARI_MODEL_HH
