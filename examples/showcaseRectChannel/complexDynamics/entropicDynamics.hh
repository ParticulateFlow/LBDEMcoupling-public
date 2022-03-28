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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef ENTROPIC_LB_DYNAMICS_HH
#define ENTROPIC_LB_DYNAMICS_HH

#include <algorithm>
#include <limits>

#include "complexDynamics/entropicDynamics.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/entropicLbTemplates.h"
#include "latticeBoltzmann/externalForceTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
struct ELBMOperations {
    static bool getAlpha(
        T &alpha, Array<T, Descriptor<T>::q> const &f, Array<T, Descriptor<T>::q> const &fNeq)
    {
        typedef Descriptor<T> L;

        T dev = 0.0;
        T t, p, inv;

        for (plint iPop = 0; iPop < L::q; ++iPop) {
            t = (T)std::fabs(fNeq[iPop] / f[iPop]);
            if (t > dev)
                dev = t;
            if (dev > 0.01)
                break;
        }

        if (dev < 1.e-6) {
            alpha = (T)2.;
            return true;
        } else if (dev < 0.01) {
            T a1 = T();
            T a2 = T();
            T a3 = T();
            T a4 = T();

            for (plint iPop = 0; iPop < L::q; ++iPop) {
                t = fNeq[iPop];
                inv = (T)1.0 / f[iPop];
                p = t * t * inv;
                t *= inv;
                a1 += p;
                p = p * t;
                a2 += p;
                p = p * t;
                a3 += p;
                p = p * t;
                a4 += p;
            }
            a1 *= 0.5;
            a2 *= -(T)1.0 / (T)6.0;
            a3 *= (T)1.0 / (T)12.0;
            a4 *= (T)-0.05;

            alpha = 2 - 4 * a2 / a1 + 16 * a2 * a2 / (a1 * a1) - 8 * a3 / a1
                    + 80 * a2 * a3 / (a1 * a1) - 80 * a2 * a2 * a2 / (a1 * a1 * a1) - 16 * a4 / a1;
            return true;
        } else {
            const T epsilon = std::numeric_limits<T>::epsilon();

            T alphaGuess = T();
            const T var = 100.0;
            const T errorMax = epsilon * var;
            T error = 1.0;
            plint count = 0;
            for (count = 0; count < 10000; ++count) {
                T entGrowth = computeEntropyGrowth(f, fNeq, alpha);
                T entGrowthDerivative = computeEntropyGrowthDerivative(f, fNeq, alpha);
                if ((error < errorMax) || (std::fabs(entGrowth) < var * epsilon)) {
                    return true;
                }
                alphaGuess = alpha - entGrowth / entGrowthDerivative;
                error = std::fabs(alpha - alphaGuess);
                alpha = alphaGuess;
            }
            return false;
        }
    }
    static T computeEntropyGrowthDerivative(
        Array<T, Descriptor<T>::q> const &f, Array<T, Descriptor<T>::q> const &fNeq, T alpha)
    {
        typedef Descriptor<T> L;

        T entropyGrowthDerivative = T();
        for (plint iPop = 0; iPop < L::q; ++iPop) {
            T tmp = f[iPop] - alpha * fNeq[iPop];
            PLB_ASSERT(tmp > T());
            entropyGrowthDerivative += fNeq[iPop] * ((T)1 + log(tmp / L::t[iPop]));
        }

        return entropyGrowthDerivative;
    }
    static T computeEntropyGrowth(
        Array<T, Descriptor<T>::q> const &f, Array<T, Descriptor<T>::q> const &fNeq, T alpha)
    {
        typedef Descriptor<T> L;

        Array<T, Descriptor<T>::q> fAlphaFneq;
        for (plint iPop = 0; iPop < L::q; ++iPop) {
            fAlphaFneq[iPop] = f[iPop] - alpha * fNeq[iPop];
        }

        return computeEntropy(f) - computeEntropy(fAlphaFneq);
    }
    static T computeEntropy(Array<T, Descriptor<T>::q> const &f)
    {
        typedef Descriptor<T> L;
        T entropy = T();
        for (plint iPop = 0; iPop < L::q; ++iPop) {
            PLB_ASSERT(f[iPop] > T());
            entropy += f[iPop] * log(f[iPop] / L::t[iPop]);
        }

        return entropy;
    }
};

//====================================================================//
//////////////////// Class ForcedEntropicDynamics //////////////////////
//====================================================================//

template <typename T, template <typename U> class Descriptor>
int ForcedEntropicDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, ForcedEntropicDynamics<T, Descriptor> >(
        "Entropic_Forced");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
ForcedEntropicDynamics<T, Descriptor>::ForcedEntropicDynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
ForcedEntropicDynamics<T, Descriptor>::ForcedEntropicDynamics(
    HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
ForcedEntropicDynamics<T, Descriptor> *ForcedEntropicDynamics<T, Descriptor>::clone() const
{
    return new ForcedEntropicDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int ForcedEntropicDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
T ForcedEntropicDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    //    T rho = Descriptor<T>::fullRho(rhoBar);
    T invRho = Descriptor<T>::invRho(rhoBar);
    Array<T, Descriptor<T>::d> u;
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        u[iD] = j[iD] * invRho;
    }
    return entropicLbTemplates<T, Descriptor>::equilibrium(iPop, rhoBar, u);
}

template <typename T, template <typename U> class Descriptor>
void ForcedEntropicDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef Descriptor<T> L;
    typedef entropicLbTemplates<T, Descriptor> eLbTempl;

    T rho;
    Array<T, Descriptor<T>::d> u;
    momentTemplates<T, Descriptor>::compute_rho_uLb(cell, rho, u);
    T uSqr = VectorTemplate<T, Descriptor>::normSqr(u);

    Array<T, Descriptor<T>::q> f, fEq, fNeq;
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        fEq[iPop] = eLbTempl::equilibrium(iPop, rho, u);
        fNeq[iPop] = cell[iPop] - fEq[iPop];
        f[iPop] = cell[iPop] + L::t[iPop];
        fEq[iPop] += L::t[iPop];
    }
    //==============================================================================//
    //============= Evaluation of alpha using a Newton Raphson algorithm ===========//
    //==============================================================================//

    T alpha = (T)2;
#ifdef PLB_DEBUG
    bool converged =
#endif
        getAlpha(alpha, f, fNeq);

    PLB_ASSERT(converged);

    T *force = cell.getExternal(forceBeginsAt);
    for (int iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
        u[iDim] += force[iDim] / (T)2.;
    }
    uSqr = VectorTemplate<T, Descriptor>::normSqr(u);
    T omegaTot = this->getOmega() / (T)2 * alpha;
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] *= (T)1 - omegaTot;
        cell[iPop] += omegaTot * eLbTempl::equilibrium(iPop, rho, u);
    }
    externalForceTemplates<T, Descriptor>::addGuoForce(cell, u, omegaTot, (T)1);

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, Descriptor<T>::rhoBar(rho), uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T ForcedEntropicDynamics<T, Descriptor>::computeEntropy(Array<T, Descriptor<T>::q> const &f)
{
    typedef Descriptor<T> L;
    T entropy = T();
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        PLB_ASSERT(f[iPop] > T());
        entropy += f[iPop] * log(f[iPop] / L::t[iPop]);
    }

    return entropy;
}

template <typename T, template <typename U> class Descriptor>
T ForcedEntropicDynamics<T, Descriptor>::computeEntropyGrowth(
    Array<T, Descriptor<T>::q> const &f, Array<T, Descriptor<T>::q> const &fNeq, T alpha)
{
    typedef Descriptor<T> L;

    Array<T, Descriptor<T>::q> fAlphaFneq;
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        fAlphaFneq[iPop] = f[iPop] - alpha * fNeq[iPop];
    }

    return computeEntropy(f) - computeEntropy(fAlphaFneq);
}

template <typename T, template <typename U> class Descriptor>
T ForcedEntropicDynamics<T, Descriptor>::computeEntropyGrowthDerivative(
    Array<T, Descriptor<T>::q> const &f, Array<T, Descriptor<T>::q> const &fNeq, T alpha)
{
    typedef Descriptor<T> L;

    T entropyGrowthDerivative = T();
    for (plint iPop = 0; iPop < L::q; ++iPop) {
        T tmp = f[iPop] - alpha * fNeq[iPop];
        PLB_ASSERT(tmp > T());
        entropyGrowthDerivative += fNeq[iPop] * ((T)1 + log(tmp / L::t[iPop]));
    }

    return entropyGrowthDerivative;
}

template <typename T, template <typename U> class Descriptor>
bool ForcedEntropicDynamics<T, Descriptor>::getAlpha(
    T &alpha, Array<T, Descriptor<T>::q> const &f, Array<T, Descriptor<T>::q> const &fNeq)
{
    const T epsilon = std::numeric_limits<T>::epsilon();

    T alphaGuess = T();
    const T var = 100.0;
    const T errorMax = epsilon * var;
    T error = 1.0;
    plint count = 0;
    for (count = 0; count < 10000; ++count) {
        T entGrowth = computeEntropyGrowth(f, fNeq, alpha);
        T entGrowthDerivative = computeEntropyGrowthDerivative(f, fNeq, alpha);
        if ((error < errorMax) || (std::fabs(entGrowth) < var * epsilon)) {
            return true;
        }
        alphaGuess = alpha - entGrowth / entGrowthDerivative;
        error = std::fabs(alpha - alphaGuess);
        alpha = alphaGuess;
    }
    return false;
}

//============================================================================//
/////////////////////////// Class EntropicDynamics ///////////////////////
//============================================================================//

template <typename T, template <typename U> class Descriptor>
int EntropicDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, EntropicDynamics<T, Descriptor> >("Entropic");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param moments_ a Moments object to know how to compute velocity moments
 */
template <typename T, template <typename U> class Descriptor>
EntropicDynamics<T, Descriptor>::EntropicDynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
EntropicDynamics<T, Descriptor>::EntropicDynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
EntropicDynamics<T, Descriptor> *EntropicDynamics<T, Descriptor>::clone() const
{
    return new EntropicDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int EntropicDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
T EntropicDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    //    T rho = Descriptor<T>::fullRho(rhoBar);
    T invRho = Descriptor<T>::invRho(rhoBar);
    Array<T, Descriptor<T>::d> u;
    for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
        u[iD] = j[iD] * invRho;
    }
    return entropicLbTemplates<T, Descriptor>::equilibrium(iPop, rhoBar, u);

    //    T invRho = Descriptor<T>::invRho(rhoBar);
    //    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j,
    //    jSqr);
}

template <typename T, template <typename U> class Descriptor>
void EntropicDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    typedef Descriptor<T> L;

    Array<T, Descriptor<T>::q> f, fEq, fNeq;
    T rhoBar;
    Array<T, Descriptor<T>::d> j;

    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    const T jSqr = VectorTemplateImpl<T, Descriptor<T>::d>::normSqr(j);

    for (plint iPop = 0; iPop < L::q; ++iPop) {
        fEq[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr);
        fNeq[iPop] = cell[iPop] - fEq[iPop];
        f[iPop] = cell[iPop] + L::t[iPop];
    }

    T alpha = (T)2;
#ifdef PLB_DEBUG
    bool converged =
#endif
        ELBMOperations<T, Descriptor>::getAlpha(alpha, f, fNeq);

    PLB_ASSERT(converged);

    T omegaTot = this->getOmega() / (T)2 * alpha;

    for (plint iPop = 0; iPop < L::q; ++iPop) {
        cell[iPop] *= (T)1 - omegaTot;
        cell[iPop] += omegaTot * fEq[iPop];
    }

    T invRho = L::invRho(rhoBar);
    T uSqr = jSqr * invRho * invRho;

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void EntropicDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    typedef Descriptor<T> L;

    Array<T, Descriptor<T>::q> f, fEq, fNeq;

    const T jSqr = VectorTemplateImpl<T, Descriptor<T>::d>::normSqr(j);

    for (plint iPop = 0; iPop < L::q; ++iPop) {
        fEq[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr);
        fNeq[iPop] = cell[iPop] - fEq[iPop];
        f[iPop] = cell[iPop] + L::t[iPop];
    }

    T alpha = (T)2;
#ifdef PLB_DEBUG
    bool converged =
#endif
        ELBMOperations<T, Descriptor>::getAlpha(alpha, f, fNeq);

    PLB_ASSERT(converged);

    T omegaTot = this->getOmega() / (T)2 * alpha;

    for (plint iPop = 0; iPop < L::q; ++iPop) {
        cell[iPop] *= (T)1 - omegaTot;
        cell[iPop] += omegaTot * fEq[iPop];
    }

    T invRho = L::invRho(rhoBar);
    T uSqr = jSqr * invRho * invRho;

    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

//============================================================================//
//////////////////////// Class VariableOmegaELBMDynamics ///////////////////////
//============================================================================//
template <typename T, template <typename U> class Descriptor>
int VariableOmegaELBMDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, VariableOmegaELBMDynamics<T, Descriptor> >(
        "VariableOmegaELBMDynamics");

template <typename T, template <typename U> class Descriptor>
VariableOmegaELBMDynamics<T, Descriptor>::VariableOmegaELBMDynamics(
    T omega0_, bool automaticPrepareCollision_) :
    VariableOmegaDynamics<T, Descriptor>(
        new BGKdynamics<T, Descriptor>(omega0_), automaticPrepareCollision_),
    omega0(omega0_)
{ }

template <typename T, template <typename U> class Descriptor>
VariableOmegaELBMDynamics<T, Descriptor>::VariableOmegaELBMDynamics(
    HierarchicUnserializer &unserializer) :
    VariableOmegaDynamics<T, Descriptor>(0, false), omega0(T())
{
    unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
VariableOmegaELBMDynamics<T, Descriptor> *VariableOmegaELBMDynamics<T, Descriptor>::clone() const
{
    return new VariableOmegaELBMDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int VariableOmegaELBMDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void VariableOmegaELBMDynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    serializer.addValue(omega0);
    VariableOmegaDynamics<T, Descriptor>::serialize(serializer);
}

template <typename T, template <typename U> class Descriptor>
void VariableOmegaELBMDynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    unserializer.readValue(omega0);
    VariableOmegaDynamics<T, Descriptor>::unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void VariableOmegaELBMDynamics<T, Descriptor>::setOmega(T omega_)
{
    omega0 = omega_;
    this->getBaseDynamics().setOmega(omega0);
}

template <typename T, template <typename U> class Descriptor>
T VariableOmegaELBMDynamics<T, Descriptor>::getOmega() const
{
    return omega0;
}

template <typename T, template <typename U> class Descriptor>
T VariableOmegaELBMDynamics<T, Descriptor>::getOmegaFromCell(Cell<T, Descriptor> const &cell) const
{
    typedef Descriptor<T> L;

    Array<T, L::q> f, fEq, fNeq;
    T rhoBar;
    Array<T, L::d> j;

    this->computeRhoBarJ(cell, rhoBar, j);
    const T jSqr = VectorTemplateImpl<T, L::d>::normSqr(j);

    for (plint iPop = 0; iPop < L::q; ++iPop) {
        fEq[iPop] = this->computeEquilibrium(iPop, rhoBar, j, jSqr);
        fNeq[iPop] = cell[iPop] - fEq[iPop];
        f[iPop] = cell[iPop] + L::t[iPop];
    }

    T alpha = (T)2;
#ifdef PLB_DEBUG
    bool converged =
#endif
        ELBMOperations<T, Descriptor>::getAlpha(alpha, f, fNeq);

    PLB_ASSERT(converged);

    return omega0 / (T)2 * alpha;
}

}  // namespace plb

#endif
