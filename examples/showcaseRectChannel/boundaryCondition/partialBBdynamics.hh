#ifndef PARTIAL_BB_DYNAMICS_HH
#define PARTIAL_BB_DYNAMICS_HH

#include <math.h>

#include <algorithm>
#include <limits>

#include "partialBBdynamics.h"

namespace plb {
/* *************** Class PartialBBdynamics *********************************************** */

template <typename T, template <typename U> class Descriptor>
int PartialBBdynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, PartialBBdynamics<T, Descriptor>>("Partial_BB");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
PartialBBdynamics<T, Descriptor>::PartialBBdynamics(T omega_) :
    IsoThermalBulkDynamics<T, Descriptor>(omega_)
{ }

template <typename T, template <typename U> class Descriptor>
PartialBBdynamics<T, Descriptor>::PartialBBdynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
PartialBBdynamics<T, Descriptor> *PartialBBdynamics<T, Descriptor>::clone() const
{
    return new PartialBBdynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int PartialBBdynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void PartialBBdynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{  // get solid fraction of the cell
    T solidFraction = *cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt);
    // and wall velocity
    Array<T, Descriptor<T>::d> wallVelocityFromExternal;
    wallVelocityFromExternal.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::jBeginsAt));

    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr =
        PSMCollision(cell, rhoBar, j, wallVelocityFromExternal, solidFraction, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void PartialBBdynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    PLB_ASSERT(true);
}

template <typename T, template <typename U> class Descriptor>
T PartialBBdynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template <typename T, template <typename U> class Descriptor>
void PartialBBdynamics<T, Descriptor>::decomposeOrder0(
    Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);

    Array<T, Descriptor<T>::q> fEq;
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(
        rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[1 + Descriptor<T>::d + iPop] = cell[iPop] - fEq[iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset + iExt] = *cell.getExternal(iExt);
    }
}

template <typename T, template <typename U> class Descriptor>
void PartialBBdynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, Descriptor<T>::q> fEq;
    dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibria(
        rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }

    int offset = 1 + Descriptor<T>::d + Descriptor<T>::q;
    for (plint iExt = 0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset + iExt];
    }
}

template <typename T, template <typename U> class Descriptor>
T PartialBBdynamics<T, Descriptor>::PSMCollision(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
    Array<T, Descriptor<T>::d> const &wallVelocity, T &solidFraction, T omega)
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    const T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    // used with Krueger method
    // T fullRho = Descriptor<T>::fullRho(rhoBar);
    // const T jSqrWall = VectorTemplate<T,Descriptor>::normSqr(fullRho*wallVelocity);

    // recover populations for solid collision
    Array<T, Descriptor<T>::q> fTmp;
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        fTmp[iPop] = cell[iPop];
    }

    // Weighting parameter for collision (fluid or solid) (see Krueger, Partially Saturated BB
    // chapter)
    T B = solidFraction;  // in Walsh model : B = ns ; solidFraction = ns*

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        plint iPopOpp = indexTemplates::opposite<Descriptor<T>>(iPop);

        // ------------- Krueger method with B = sF------------
        // // Fluid collision
        // // f_i - f_i*(1-B)/tau  1st part of eq 11.10 Krueger
        // cell[iPop] *= (T)1-omega*(T(1)-B);
        // // (1-B)*feq/tau    2nd part of eq 11.10 Krueger
        // cell[iPop] += omega*(T(1)-B)*dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop,
        // rhoBar, invRho, j, jSqr);

        // // Solid collision
        // cell[iPop] += B*(fTmp[iPopOpp] -
        // dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPopOpp,rhoBar, invRho, j, jSqr));
        // cell[iPop] -= B*(fTmp[iPop] -
        // dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop,rhoBar, invRho,
        // fullRho*wallVelocity, jSqrWall));

        // ----------- Walsh method "a new pbb model" eq 8 ------------------------
        // Notation : code = Krueger = Walsh
        // fTmp[iPop] = f_i(x,t) = f^in_a(x,t)
        // cell[iPop] = f_i (x+c_i*dt, t+dt)
        // bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr) = f_i eq(x,t) = feq(x,t)
        // fTmp[iPop] + omega*(-fTmp[iPop]+bgk_ma2_equi)  = | f_i(x,t) + (-f_i(x,t)+feq_i(x,t))/tau
        // = f_i*(x,t) | = f^c_a (x,t) : f just after fluid collision

        cell[iPop] = (T(1.) - B)
                         * (cell[iPop] * (1. - omega)
                            + omega
                                  * dynamicsTemplates<T, Descriptor>::bgk_ma2_equilibrium(
                                      iPop, rhoBar, invRho, j, jSqr))
                     + B * fTmp[iPopOpp];
    }
    return invRho * invRho * jSqr;  // return uSqr
}

}  // namespace plb

#endif  // PARTIAL_BB_DYNAMICS_HH
