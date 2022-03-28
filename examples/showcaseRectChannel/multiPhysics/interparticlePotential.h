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

#ifndef INTERPARTICLE_POTENTIAL_H
#define INTERPARTICLE_POTENTIAL_H

#include <cmath>

#include "core/globalDefs.h"

namespace plb {

/// Interaction potentials which depend only on the density rho
namespace interparticlePotential {

/// Interface for the interaction potential ("psi function") for single-component multi-phase
/// Shan/Chen
template <typename T>
struct PsiFunction {
    virtual ~PsiFunction();
    virtual T compute(T rho) const = 0;
    virtual PsiFunction<T> *clone() const = 0;
};

template <typename T>
class PsiIsRho : public PsiFunction<T> {
    virtual T compute(T rho) const;
    virtual PsiIsRho<T> *clone() const;
};

/// Shan, Chen 1993: /f[ \Psi = \rho_0(1-\exp(-\rho/\rho_0)) /f]
template <typename T>
class PsiShanChen93 : public PsiFunction<T> {
public:
    PsiShanChen93(T rho_0_);
    virtual T compute(T rho) const;
    virtual PsiShanChen93<T> *clone() const;

private:
    T rho_0;
};

/// Shan, Chen 1994: /f[ \Psi = \Psi_0\exp(-\rho_0/\rho) /f]
template <typename T>
class PsiShanChen94 : public PsiFunction<T> {
public:
    PsiShanChen94(T psi_0_, T rho_0_);
    virtual T compute(T rho) const;
    virtual PsiShanChen94<T> *clone() const;

private:
    T psi_0, rho_0;
};

/// Qian e.a. 1995: \f[ \Psi = g\rho_0^2\rho^2/(2(\rho_0+\rho)^2) \f]
template <typename T>
class PsiQian95 : public PsiFunction<T> {
public:
    PsiQian95(T rho_0_, T g_);
    virtual T compute(T rho) const;
    virtual PsiQian95<T> *clone() const;

private:
    T rho_0, rho_0_sqr, g;
};

}  // namespace interparticlePotential

}  // namespace plb

#endif  // INTERPARTICLE_POTENTIAL_H
