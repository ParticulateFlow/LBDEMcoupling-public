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

#ifndef TRT_DYNAMICS_HH
#define TRT_DYNAMICS_HH

#include <algorithm>
#include <limits>

#include "complexDynamics/trtDynamics.h"
#include "core/latticeStatistics.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

/* *************** Class BaseTRTdynamics *********************************************** */

/** \param omegaPlus_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
BaseTRTdynamics<T, Descriptor>::BaseTRTdynamics(T omegaPlus_, T omegaMinus_, bool constant_magic) :
    IsoThermalBulkDynamics<T, Descriptor>(omegaPlus_),
    keep_magic_constant_when_setting_omega(constant_magic)
{
    if (!util::fpequal(omegaMinus_, (T)0.0))
        omegaMinus = omegaMinus_;
    else
        omegaMinus = (T)8 * ((T)2 - omegaPlus_) / ((T)8 - omegaPlus_);
}

// construct the class and initialize it using the serializer
template <typename T, template <typename U> class Descriptor>
BaseTRTdynamics<T, Descriptor>::BaseTRTdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(
        T()),  // dummy init, true variables values will be set by the unserializer
    omegaMinus(T()),
    keep_magic_constant_when_setting_omega(bool())
{
    this->unserialize(unserializer);
}

/**
 * Here, we need to add the additional members of the TRT class to the serializer.
 * @tparam T
 * @tparam Descriptor
 * @param serializer
 */
template <typename T, template <typename U> class Descriptor>
void BaseTRTdynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    serializer.addValue(omegaMinus);
    serializer.addValue(keep_magic_constant_when_setting_omega);
    IsoThermalBulkDynamics<T, Descriptor>::serialize(serializer);
}

/**
 * Here, we need to read the value of keep_magic_constant_when_setting_omega BEFORE unserializing
 * the object, because Dynamics<T, Descriptor>::unserialize calls
 * this->setOmega(unserializer.readValue<T>()), which reads keep_magic_constant_when_setting_omega.
 * If keep_magic_constant_when_setting_omega is not initizlized, this leads to undefined behavior.
 * @tparam T
 * @tparam Descriptor
 * @param unserializer
 */
template <typename T, template <typename U> class Descriptor>
void BaseTRTdynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    omegaMinus = unserializer.readValue<T>();
    keep_magic_constant_when_setting_omega = unserializer.readValue<bool>();
    IsoThermalBulkDynamics<T, Descriptor>::unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
T BaseTRTdynamics<T, Descriptor>::getOmegaMinus() const
{
    return this->omegaMinus;
}

template <typename T, template <typename U> class Descriptor>
void BaseTRTdynamics<T, Descriptor>::setOmegaMinus(T omegaMinus_)
{
    this->omegaMinus = omegaMinus_;
}
template <typename T, template <typename U> class Descriptor>
T BaseTRTdynamics<T, Descriptor>::getMagicParam() const
{
    return (1. / this->getOmega() - 0.5) * (1. / omegaMinus - 0.5);
}

template <typename T, template <typename U> class Descriptor>
void BaseTRTdynamics<T, Descriptor>::setMagicParam(T magic_)
{
    auto omegaPlus = this->getOmega();
    this->omegaMinus = (4. - 2. * omegaPlus) / (2. - omegaPlus + 4. * magic_ * omegaPlus);
}

template <typename T, template <typename U> class Descriptor>
T BaseTRTdynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    if (whichParameter == dynamicParams::omega_plus) {
        return this->getOmega();
    } else if (whichParameter == dynamicParams::omega_minus) {
        return this->getOmegaMinus();
    } else if (whichParameter == dynamicParams::magicParameter) {
        return this->getMagicParam();
    } else {
        pcout << "Can't get the parameter\n";
        abort();
    }
    return 0.;
}

template <typename T, template <typename U> class Descriptor>
void BaseTRTdynamics<T, Descriptor>::setParameter(plint whichParameter, T value)
{
    if (whichParameter == dynamicParams::omega_plus) {
        setOmega(value);
    } else if (whichParameter == dynamicParams::omega_minus) {
        setOmegaMinus(value);
    } else if (whichParameter == dynamicParams::magicParameter) {
        setMagicParam(value);
    } else {
        pcout << "Can't set the parameter\n";
        abort();
    }
}

template <typename T, template <typename U> class Descriptor>
void BaseTRTdynamics<T, Descriptor>::setOmega(T omega_)
{
    T magic = getMagicParam();
    BasicBulkDynamics<T, Descriptor>::setOmega(omega_);
    if (keep_magic_constant_when_setting_omega)
        setMagicParam(magic);
}

/* *************** Class TRTdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int TRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, TRTdynamics<T, Descriptor>>("TRT");

template <typename T, template <typename U> class Descriptor>
TRTdynamics<T, Descriptor> *TRTdynamics<T, Descriptor>::clone() const
{
    return new TRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int TRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void TRTdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    const T omegaPlus = this->getOmega();
    const T omegaMinus = this->getOmegaMinus();

    Array<T, Descriptor<T>::q> eq;
    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    Array<T, Descriptor<T>::q / 2 + 1> eq_plus, eq_minus, f_plus, f_minus;

    Array<T, Descriptor<T>::d> j;
    T rhoBar;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = normSqr(j);
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, eq);

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        eq_plus[i] = 0.5 * (eq[i] + eq[i + Descriptor<T>::q / 2]);
        eq_minus[i] = 0.5 * (eq[i] - eq[i + Descriptor<T>::q / 2]);
        f_plus[i] = 0.5 * (cell[i] + cell[i + Descriptor<T>::q / 2]);
        f_minus[i] = 0.5 * (cell[i] - cell[i + Descriptor<T>::q / 2]);
    }

    cell[0] += -omegaPlus * cell[0] + omegaPlus * eq[0];

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        cell[i] += -omegaPlus * (f_plus[i] - eq_plus[i]) - omegaMinus * (f_minus[i] - eq_minus[i]);
        cell[i + Descriptor<T>::q / 2] +=
            -omegaPlus * (f_plus[i] - eq_plus[i]) + omegaMinus * (f_minus[i] - eq_minus[i]);
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
void TRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    const T omegaPlus = this->getOmega();
    const T omegaMinus = this->getOmegaMinus();

    Array<T, Descriptor<T>::q> eq;
    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    Array<T, Descriptor<T>::q / 2 + 1> eq_plus, eq_minus, f_plus, f_minus;

    T jSqr = normSqr(j);
    T invRho = Descriptor<T>::invRho(rhoBar);
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, eq);

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        eq_plus[i] = 0.5 * (eq[i] + eq[i + Descriptor<T>::q / 2]);
        eq_minus[i] = 0.5 * (eq[i] - eq[i + Descriptor<T>::q / 2]);
        f_plus[i] = 0.5 * (cell[i] + cell[i + Descriptor<T>::q / 2]);
        f_minus[i] = 0.5 * (cell[i] - cell[i + Descriptor<T>::q / 2]);
    }

    cell[0] += -omegaPlus * cell[0] + omegaPlus * eq[0];

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        cell[i] += -omegaPlus * (f_plus[i] - eq_plus[i]) - omegaMinus * (f_minus[i] - eq_minus[i]);
        cell[i + Descriptor<T>::q / 2] +=
            -omegaPlus * (f_plus[i] - eq_plus[i]) + omegaMinus * (f_minus[i] - eq_minus[i]);
    }

    if (cell.takesStatistics()) {
        gatherStatistics(
            statistics, rhoBar,
            jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar));
    }
}

template <typename T, template <typename U> class Descriptor>
T TRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/* *************** Class Ma1TRTdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int Ma1TRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, TRTdynamics<T, Descriptor>>("Ma1TRT");

template <typename T, template <typename U> class Descriptor>
Ma1TRTdynamics<T, Descriptor> *Ma1TRTdynamics<T, Descriptor>::clone() const
{
    return new Ma1TRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int Ma1TRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void Ma1TRTdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    const T omegaPlus = this->getOmega();
    const T omegaMinus = this->getOmegaMinus();

    Array<T, Descriptor<T>::q> eq;
    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    Array<T, Descriptor<T>::q / 2 + 1> eq_plus, eq_minus, f_plus, f_minus;

    Array<T, Descriptor<T>::d> j;
    T rhoBar;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = normSqr(j);
    T invRho = Descriptor<T>::invRho(rhoBar);
    for (int i = 0; i < Descriptor<T>::q; ++i) {
        eq[i] = dynamicsTemplatesImpl<T, Descriptor<T>>::bgk_ma1_equilibrium(i, rhoBar, j);
    }

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        eq_plus[i] = 0.5 * (eq[i] + eq[i + Descriptor<T>::q / 2]);
        eq_minus[i] = 0.5 * (eq[i] - eq[i + Descriptor<T>::q / 2]);
        f_plus[i] = 0.5 * (cell[i] + cell[i + Descriptor<T>::q / 2]);
        f_minus[i] = 0.5 * (cell[i] - cell[i + Descriptor<T>::q / 2]);
    }

    cell[0] += -omegaPlus * cell[0] + omegaPlus * eq[0];

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        cell[i] += -omegaPlus * (f_plus[i] - eq_plus[i]) - omegaMinus * (f_minus[i] - eq_minus[i]);
        cell[i + Descriptor<T>::q / 2] +=
            -omegaPlus * (f_plus[i] - eq_plus[i]) + omegaMinus * (f_minus[i] - eq_minus[i]);
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr * invRho * invRho);
    }
}

template <typename T, template <typename U> class Descriptor>
void Ma1TRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    const T omegaPlus = this->getOmega();
    const T omegaMinus = this->getOmegaMinus();

    Array<T, Descriptor<T>::q> eq;
    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    Array<T, Descriptor<T>::q / 2 + 1> eq_plus, eq_minus, f_plus, f_minus;

    T jSqr = normSqr(j);

    for (int i = 0; i < Descriptor<T>::q; ++i) {
        eq[i] = dynamicsTemplatesImpl<T, Descriptor<T>>::bgk_ma1_equilibrium(i, rhoBar, j);
    }

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        eq_plus[i] = 0.5 * (eq[i] + eq[i + Descriptor<T>::q / 2]);
        eq_minus[i] = 0.5 * (eq[i] - eq[i + Descriptor<T>::q / 2]);
        f_plus[i] = 0.5 * (cell[i] + cell[i + Descriptor<T>::q / 2]);
        f_minus[i] = 0.5 * (cell[i] - cell[i + Descriptor<T>::q / 2]);
    }

    cell[0] += -omegaPlus * cell[0] + omegaPlus * eq[0];

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        cell[i] += -omegaPlus * (f_plus[i] - eq_plus[i]) - omegaMinus * (f_minus[i] - eq_minus[i]);
        cell[i + Descriptor<T>::q / 2] +=
            -omegaPlus * (f_plus[i] - eq_plus[i]) + omegaMinus * (f_minus[i] - eq_minus[i]);
    }

    if (cell.takesStatistics()) {
        gatherStatistics(
            statistics, rhoBar,
            jSqr * Descriptor<T>::invRho(rhoBar) * Descriptor<T>::invRho(rhoBar));
    }
}

template <typename T, template <typename U> class Descriptor>
T Ma1TRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    return dynamicsTemplatesImpl<T, Descriptor<T>>::bgk_ma1_equilibrium(iPop, rhoBar, j);
}

/* *************** Class IncTRTdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int IncTRTdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, IncTRTdynamics<T, Descriptor>>("IncTRT");

template <typename T, template <typename U> class Descriptor>
IncTRTdynamics<T, Descriptor> *IncTRTdynamics<T, Descriptor>::clone() const
{
    return new IncTRTdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int IncTRTdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void IncTRTdynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    const T omegaPlus = this->getOmega();
    const T omegaMinus = this->getOmegaMinus();

    Array<T, Descriptor<T>::q> eq;
    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    Array<T, Descriptor<T>::q / 2 + 1> eq_plus, eq_minus, f_plus, f_minus;

    Array<T, 3> j;
    T rhoBar;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = normSqr(j);
    T invRho0 = 1.;
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(rhoBar, invRho0, j, jSqr, eq);

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        eq_plus[i] = 0.5 * (eq[i] + eq[i + Descriptor<T>::q / 2]);
        eq_minus[i] = 0.5 * (eq[i] - eq[i + Descriptor<T>::q / 2]);
        f_plus[i] = 0.5 * (cell[i] + cell[i + Descriptor<T>::q / 2]);
        f_minus[i] = 0.5 * (cell[i] - cell[i + Descriptor<T>::q / 2]);
    }

    cell[0] += -omegaPlus * cell[0] + omegaPlus * eq[0];

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        cell[i] += -omegaPlus * (f_plus[i] - eq_plus[i]) - omegaMinus * (f_minus[i] - eq_minus[i]);
        cell[i + Descriptor<T>::q / 2] +=
            -omegaPlus * (f_plus[i] - eq_plus[i]) + omegaMinus * (f_minus[i] - eq_minus[i]);
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void IncTRTdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    const T omegaPlus = this->getOmega();
    const T omegaMinus = this->getOmegaMinus();

    Array<T, Descriptor<T>::q> eq;
    // In the following, we number the plus/minus variables from 1 to (Q-1)/2.
    // So we allocate the index-zero memory location, and waste some memory
    // for convenience.
    Array<T, Descriptor<T>::q / 2 + 1> eq_plus, eq_minus, f_plus, f_minus;

    T jSqr = normSqr(j);
    T invRho0 = 1.;
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(rhoBar, invRho0, j, jSqr, eq);

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        eq_plus[i] = 0.5 * (eq[i] + eq[i + Descriptor<T>::q / 2]);
        eq_minus[i] = 0.5 * (eq[i] - eq[i + Descriptor<T>::q / 2]);
        f_plus[i] = 0.5 * (cell[i] + cell[i + Descriptor<T>::q / 2]);
        f_minus[i] = 0.5 * (cell[i] - cell[i + Descriptor<T>::q / 2]);
    }

    cell[0] += -omegaPlus * cell[0] + omegaPlus * eq[0];

    for (plint i = 1; i <= Descriptor<T>::q / 2; ++i) {
        cell[i] += -omegaPlus * (f_plus[i] - eq_plus[i]) - omegaMinus * (f_minus[i] - eq_minus[i]);
        cell[i + Descriptor<T>::q / 2] +=
            -omegaPlus * (f_plus[i] - eq_plus[i]) + omegaMinus * (f_minus[i] - eq_minus[i]);
    }

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, jSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T IncTRTdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, (T)1, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
bool IncTRTdynamics<T, Descriptor>::velIsJ() const
{
    return true;
}

template <typename T, template <typename U> class Descriptor>
void IncTRTdynamics<T, Descriptor>::computeVelocity(
    Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &u) const
{
    T dummyRhoBar;
    this->computeRhoBarJ(cell, dummyRhoBar, u);
}

}  // namespace plb

#endif  // TRT_DYNAMICS_HH
