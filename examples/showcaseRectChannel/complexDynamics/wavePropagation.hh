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
#ifndef WAVE_PROPAGATION_HH
#define WAVE_PROPAGATION_HH

#include "basicDynamics/isoThermalDynamics.h"
#include "core/cell.h"
#include "core/dynamicsIdentifiers.h"
#include "core/latticeStatistics.h"

namespace plb {

/* *************** Class WaveDynamics ************************************ */

template <typename T, template <typename U> class Descriptor>
int WaveDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, WaveDynamics<T, Descriptor> >("Wave");

/** \param vs2_ speed of sound
 *  \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template <typename T, template <typename U> class Descriptor>
WaveDynamics<T, Descriptor>::WaveDynamics(T vs2_) :
    IsoThermalBulkDynamics<T, Descriptor>((T)2.0), vs2(vs2_)
{ }

template <typename T, template <typename U> class Descriptor>
WaveDynamics<T, Descriptor>::WaveDynamics(HierarchicUnserializer &unserializer) :
    IsoThermalBulkDynamics<T, Descriptor>((T)2.0), vs2(T())
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void WaveDynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    serializer.addValue(vs2);
}

template <typename T, template <typename U> class Descriptor>
void WaveDynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    vs2 = unserializer.readValue<T>();
}

template <typename T, template <typename U> class Descriptor>
WaveDynamics<T, Descriptor> *WaveDynamics<T, Descriptor>::clone() const
{
    return new WaveDynamics<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
int WaveDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void WaveDynamics<T, Descriptor>::collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T uSqr = waveCollision(cell, rhoBar, j, vs2);
    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
void WaveDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &stat)
{
    T uSqr = waveCollision(cell, rhoBar, j, vs2);
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template <typename T, template <typename U> class Descriptor>
T WaveDynamics<T, Descriptor>::computeEquilibrium(
    plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return waveEquilibrium(iPop, rhoBar, invRho, j, jSqr, vs2);
}

template <typename T, template <typename U> class Descriptor>
T WaveDynamics<T, Descriptor>::getParameter(plint whichParameter) const
{
    switch (whichParameter) {
    case dynamicParams::sqrSpeedOfSound:
        return this->getVs2();
    };
    return 0.;
}

template <typename T, template <typename U> class Descriptor>
void WaveDynamics<T, Descriptor>::setParameter(plint whichParameter, T value)
{
    switch (whichParameter) {
    case dynamicParams::sqrSpeedOfSound:
        setVs2(value);
    };
}

template <typename T, template <typename U> class Descriptor>
T WaveDynamics<T, Descriptor>::getVs2() const
{
    return vs2;
}

template <typename T, template <typename U> class Descriptor>
void WaveDynamics<T, Descriptor>::setVs2(T vs2_)
{
    vs2 = vs2_;
}

template <typename T, template <typename U> class Descriptor>
T WaveDynamics<T, Descriptor>::waveCollision(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T vs2)
{
    const T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T invRho = Descriptor<T>::invRho(rhoBar);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] *= (T)-1.0;
        cell[iPop] += (T)2.0 * waveEquilibrium(iPop, rhoBar, invRho, j, jSqr, vs2);
    }
    return invRho * invRho * jSqr;
}

template <typename T, template <typename U> class Descriptor>
T WaveDynamics<T, Descriptor>::waveEquilibrium(
    plint iPop, T rhoBar, T invRho, Array<T, Descriptor<T>::d> const &j, T jSqr, T vs2)
{
    T kappa = vs2 - Descriptor<T>::cs2;
    if (iPop == 0) {
        return Descriptor<T>::invCs2
               * (kappa * (Descriptor<T>::t[0] - (T)1)
                  + rhoBar * (Descriptor<T>::t[0] * vs2 - kappa));
    } else {
        T c_j = T();
        for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
            c_j += Descriptor<T>::c[iPop][iD] * j[iD];
        }
        return Descriptor<T>::invCs2 * Descriptor<T>::t[iPop] * (kappa + rhoBar * vs2 + c_j);
    }
}

// ============== Wave Absorption dynamics ====================== //

template <typename T, template <typename U> class Descriptor>
int WaveAbsorptionDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, WaveAbsorptionDynamics<T, Descriptor> >(
        "Wave_Absorption");

template <typename T, template <typename U> class Descriptor>
WaveAbsorptionDynamics<T, Descriptor>::WaveAbsorptionDynamics(
    Dynamics<T, Descriptor> *baseDynamics_) :
    CompositeDynamics<T, Descriptor>(
        baseDynamics_, false)  // false is for automaticPrepareCollision.
{ }

template <typename T, template <typename U> class Descriptor>
WaveAbsorptionDynamics<T, Descriptor>::WaveAbsorptionDynamics(
    HierarchicUnserializer &unserializer) :
    CompositeDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionDynamics<T, Descriptor>::collide(
    Cell<T, Descriptor> &cell, BlockStatistics &statistics)
{
    static const T epsilon = 1.e4 * std::numeric_limits<T>::epsilon();
    int sigmaPos = Descriptor<T>::ExternalField::sigmaBeginsAt;
    int rhoBarPos = Descriptor<T>::ExternalField::rhoBarBeginsAt;
    int uPos = Descriptor<T>::ExternalField::uBeginsAt;
    T sigma = *cell.getExternal(sigmaPos);
    T rhoBarF = *cell.getExternal(rhoBarPos);
    Array<T, Descriptor<T>::d> uF;
    uF.from_cArray(cell.getExternal(uPos));

    if (sigma < epsilon) {
        this->getBaseDynamics().collide(cell, statistics);
    } else {
        T rhoBar;
        Array<T, Descriptor<T>::d> j;
        this->getBaseDynamics().computeRhoBarJ(cell, rhoBar, j);
        T jSqr = normSqr(j);
        Array<T, Descriptor<T>::q> fEq;
        this->getBaseDynamics().computeEquilibria(fEq, rhoBar, j, jSqr);

        Array<T, Descriptor<T>::d> jF = Descriptor<T>::fullRho(rhoBarF) * uF;
        T jFsqr = normSqr(jF);
        Array<T, Descriptor<T>::q> fEqF;
        this->getBaseDynamics().computeEquilibria(fEqF, rhoBarF, jF, jFsqr);

        this->getBaseDynamics().collide(cell, statistics);

        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            cell[iPop] -= sigma * (fEq[iPop] - fEqF[iPop]);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionDynamics<T, Descriptor>::collideExternal(
    Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
    BlockStatistics &statistics)
{
    static const T epsilon = 1.e4 * std::numeric_limits<T>::epsilon();
    int sigmaPos = Descriptor<T>::ExternalField::sigmaBeginsAt;
    int rhoBarPos = Descriptor<T>::ExternalField::rhoBarBeginsAt;
    int uPos = Descriptor<T>::ExternalField::uBeginsAt;
    T sigma = *cell.getExternal(sigmaPos);
    T rhoBarF = *cell.getExternal(rhoBarPos);
    Array<T, Descriptor<T>::d> uF;
    uF.from_cArray(cell.getExternal(uPos));

    if (sigma < epsilon) {
        this->getBaseDynamics().collideExternal(cell, rhoBar, j, thetaBar, statistics);
    } else {
        T jSqr = normSqr(j);
        Array<T, Descriptor<T>::q> fEq;
        this->getBaseDynamics().computeEquilibria(fEq, rhoBar, j, jSqr);

        Array<T, Descriptor<T>::d> jF = Descriptor<T>::fullRho(rhoBarF) * uF;
        T jFsqr = normSqr(jF);
        Array<T, Descriptor<T>::q> fEqF;
        this->getBaseDynamics().computeEquilibria(fEqF, rhoBarF, jF, jFsqr);

        this->getBaseDynamics().collideExternal(cell, rhoBar, j, thetaBar, statistics);

        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            cell[iPop] -= sigma * (fEq[iPop] - fEqF[iPop]);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionDynamics<T, Descriptor>::recompose(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData, plint order) const
{
    PLB_PRECONDITION(
        (plint)rawData.size() == this->getBaseDynamics().numDecomposedVariables(order));

    if (order == 0) {
        recomposeOrder0(cell, rawData);
    } else {
        recomposeOrder1(cell, rawData);
    }
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionDynamics<T, Descriptor>::recomposeOrder0(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar = rawData[0];
    Array<T, Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);

    Array<T, Descriptor<T>::q> fEq;
    this->getBaseDynamics().computeEquilibria(fEq, rhoBar, j, jSqr);

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1 + Descriptor<T>::d + iPop];
    }
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionDynamics<T, Descriptor>::recomposeOrder1(
    Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const
{
    T rhoBar;
    Array<T, Descriptor<T>::d> j;
    Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;

    rhoBar = rawData[0];
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    PiNeq.from_cArray(&rawData[1 + Descriptor<T>::d]);

    this->getBaseDynamics().regularize(cell, rhoBar, j, jSqr, PiNeq);
}

template <typename T, template <typename U> class Descriptor>
void WaveAbsorptionDynamics<T, Descriptor>::prepareCollision(Cell<T, Descriptor> &cell)
{ }

template <typename T, template <typename U> class Descriptor>
int WaveAbsorptionDynamics<T, Descriptor>::getId() const
{
    return id;
}

// Implementation of a specific "sigma" function for WaveAbsorptionDynamics.

template <typename T>
WaveAbsorptionSigmaFunction3D<T>::WaveAbsorptionSigmaFunction3D(
    Box3D domain_, Array<plint, 6> const &numCells_, T omega_) :
    domain(domain_), numCells(numCells_), xi((T)4 / omega_ - (T)1.0e-3)
{ }

template <typename T>
T WaveAbsorptionSigmaFunction3D<T>::operator()(plint iX, plint iY, plint iZ) const
{
    std::vector<plint> distances(6, 0);

    numCells[0] ? addDistance(domain.x0 + numCells[0], iX, distances, 0) : (void)0;
    numCells[1] ? addDistance(iX, domain.x1 - numCells[1], distances, 1) : (void)0;
    numCells[2] ? addDistance(domain.y0 + numCells[2], iY, distances, 2) : (void)0;
    numCells[3] ? addDistance(iY, domain.y1 - numCells[3], distances, 3) : (void)0;
    numCells[4] ? addDistance(domain.z0 + numCells[4], iZ, distances, 4) : (void)0;
    numCells[5] ? addDistance(iZ, domain.z1 - numCells[5], distances, 5) : (void)0;

    plint distance = 0;
    plint ind = -1;
    for (pluint i = 0; i < distances.size(); ++i) {
        if (distances[i] > distance) {
            distance = distances[i];
            ind = i;
        }
    }

    if (distance == 0) {
        return (T());
    } else {
        return (xi * sigma(T(), (T)numCells[ind], (T)distance));
    }
}

template <typename T>
void WaveAbsorptionSigmaFunction3D<T>::addDistance(
    plint from, plint pos, std::vector<plint> &distances, plint i) const
{
    plint dist = from - pos;
    if (dist > 0) {
        distances[i] = dist;
    }
}

template <typename T>
T WaveAbsorptionSigmaFunction3D<T>::sigma(T x0, T x1, T x) const
{
    return ((T)3125 * (x1 - x) * std::pow(x - x0, (T)4)) / ((T)256 * std::pow(x1 - x0, (T)5));
}

// Implementation of a 2D specific "sigma" function for WaveAbsorptionDynamics.

template <typename T>
WaveAbsorptionSigmaFunction2D<T>::WaveAbsorptionSigmaFunction2D(
    Box2D domain_, Array<plint, 4> const &numCells_, T omega_) :
    domain(domain_),
    numCells(numCells_),
    xi((T)4 / omega_ - (T)1.0e-3)  // using (xi < 4/omega) or (xi < 4 * tau) to maintain the
                                   // stability xi = 4/s - epsilon, where epsilon = 0.001
{ }

template <typename T>
T WaveAbsorptionSigmaFunction2D<T>::operator()(plint iX, plint iY) const
{
    std::vector<plint> distances(4, 0);

    // if numCells[X] != 0  --> addDistance, else --> (void) 0;
    numCells[0] ? addDistance(domain.x0 + numCells[0], iX, distances, 0) : (void)0;
    numCells[1] ? addDistance(iX, domain.x1 - numCells[1], distances, 1) : (void)0;
    numCells[2] ? addDistance(domain.y0 + numCells[2], iY, distances, 2) : (void)0;
    numCells[3] ? addDistance(iY, domain.y1 - numCells[3], distances, 3) : (void)0;

    // traverse distances[i] to find the largest distance and use it to calculate the sigma
    plint distance = 0;
    plint ind = -1;
    for (pluint i = 0; i < distances.size(); ++i) {
        if (distances[i] > distance) {
            distance = distances[i];
            ind = i;
        }
    }

    if (distance == 0) {
        return (T());
    } else {
        // Kam's form
        // return 0.3 * pow((T)distance/(T)numCells[ind],2);
        // Xu and Sagaut's form
        return (xi * sigma(T(), (T)numCells[ind], (T)distance));
    }
}

template <typename T>
void WaveAbsorptionSigmaFunction2D<T>::addDistance(
    plint from, plint pos, std::vector<plint> &distances, plint i) const
{
    plint dist = from - pos;
    if (dist > 0) {
        distances[i] = dist;
    }
}

template <typename T>
T WaveAbsorptionSigmaFunction2D<T>::sigma(T x0, T x1, T x) const
{
    return ((T)3125 * (x1 - x) * std::pow(x - x0, (T)4)) / ((T)256 * std::pow(x1 - x0, (T)5));

    /*            3125(L-x)(x-x0)^4
     * sigma_x = -------------------
     *               256(L-x0)^5
     * , where x0 = 0, L or x1 = numCells[ind] and x = distance.
     */
}

}  // namespace plb

#endif  // WAVE_PROPAGATION_HH
