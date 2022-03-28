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

#ifndef INTERPARTICLE_POTENTIAL_HH
#define INTERPARTICLE_POTENTIAL_HH

#include "multiPhysics/interparticlePotential.h"

namespace plb {

namespace interparticlePotential {

template <typename T>
PsiFunction<T>::~PsiFunction()
{ }

template <typename T>
T PsiIsRho<T>::compute(T rho) const
{
    return rho;
}
template <typename T>
PsiIsRho<T> *PsiIsRho<T>::clone() const
{
    return new PsiIsRho<T>(*this);
}

template <typename T>
PsiShanChen93<T>::PsiShanChen93(T rho_0_) : rho_0(rho_0_)
{ }

template <typename T>
T PsiShanChen93<T>::compute(T rho) const
{
    return rho_0 * ((T)1 - std::exp(-rho / rho_0));
}

template <typename T>
PsiShanChen93<T> *PsiShanChen93<T>::clone() const
{
    return new PsiShanChen93<T>(*this);
}

template <typename T>
PsiShanChen94<T>::PsiShanChen94(T psi_0_, T rho_0_) : psi_0(psi_0_), rho_0(rho_0_)
{ }

template <typename T>
T PsiShanChen94<T>::compute(T rho) const
{
    return psi_0 * std::exp(-rho_0 / rho);
}

template <typename T>
PsiShanChen94<T> *PsiShanChen94<T>::clone() const
{
    return new PsiShanChen94<T>(*this);
}

template <typename T>
PsiQian95<T>::PsiQian95(T rho_0_, T g_) : rho_0(rho_0_), rho_0_sqr(rho_0 * rho_0), g(g_)
{ }

template <typename T>
T PsiQian95<T>::compute(T rho) const
{
    T rho_sqr = rho * rho;
    return g * rho_0_sqr * rho_sqr / ((T)2 * (rho_0_sqr + rho_sqr + (T)2 * rho_0 * rho));
}
template <typename T>
PsiQian95<T> *PsiQian95<T>::clone() const
{
    return new PsiQian95<T>(*this);
}

}  // namespace interparticlePotential

}  // namespace plb

#endif  // INTERPARTICLE_POTENTIAL_HH
