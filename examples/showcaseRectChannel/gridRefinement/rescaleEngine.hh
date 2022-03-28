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
 * Coupling between grids of different refinement level -- generic implementation.
 */
#ifndef GRID_REFINEMENT_HELPERS_HH
#define GRID_REFINEMENT_HELPERS_HH

#include "complexDynamics/smagorinskyDynamics.hh"
#include "gridRefinement/rescaleEngine.h"

namespace plb {

/* ====================== Rescaler base class ===================*/
template <typename T, template <typename U> class Descriptor>
Rescaler<T, Descriptor>::Rescaler()
{ }

template <typename T, template <typename U> class Descriptor>
Rescaler<T, Descriptor>::~Rescaler()
{ }

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::q> Rescaler<T, Descriptor>::computeRescaledRelFreq(
    const Array<T, Descriptor<T>::q> &relFreq, T xDt) const
{
    Array<T, Descriptor<T>::q> resRelFreq;
    //     T oneOverDt = (T)1/xDt;
    //     T twoXdt = (T)2*xDt;
    //     T xDt_1 = xDt-(T)1;

    for (pluint iA = 0; iA < Descriptor<T>::q; ++iA) {
        // resRelFreq[iA] = (T)1/(((T)1/relFreq[iA]-(T)0.5)/xDt+(T)0.5);
        resRelFreq[iA] = xDt * relFreq[iA] / ((T)1 - 0.5 * relFreq[iA] * ((T)1 - xDt));
    }

    return resRelFreq;
}

template <typename T, template <typename U> class Descriptor>
void Rescaler<T, Descriptor>::decomposeAndRescale(
    Cell<T, Descriptor> const &cell, T xDt, plint order, std::vector<T> &decompAndRescaled) const
{
    PLB_ASSERT(order == 0 || order == 1);
    cell.getDynamics().decompose(cell, decompAndRescaled, order);
    rescale(cell.getDynamics(), xDt, decompAndRescaled);
}

/* ====================== Convective no force rescaler class ===================*/
template <typename T, template <typename U> class Descriptor>
ConvectiveNoForceRescaler<T, Descriptor>::ConvectiveNoForceRescaler()
{ }

template <typename T, template <typename U> class Descriptor>
void ConvectiveNoForceRescaler<T, Descriptor>::rescale(
    const Dynamics<T, Descriptor> &dyn, T xDt, std::vector<T> &rawData) const
{
    Array<T, Descriptor<T>::q> relFreq = dyn.getRelaxationFrequencies();

    // rescale value of the relaxation time
    Array<T, Descriptor<T>::q> resRelFreq = this->computeRescaledRelFreq(relFreq, xDt);
    // rawData[0] = rho => no rescale
    // rawData[1]-rawData[1+d] => velocity => no rescale in the convective rescale
    // only the non equilibrium part is rescaled
    T prefactor = relFreq[0] / resRelFreq[0] * xDt;
    for (plint iPop = 0; iPop
                         < ((rawData.size() - Descriptor<T>::ExternalField::numScalars)
                            - ((T)1 + Descriptor<T>::d));
         ++iPop)
    {
        rawData[1 + Descriptor<T>::d + iPop] *= prefactor;
    }
}

template <typename T, template <typename U> class Descriptor>
ConvectiveNoForceRescaler<T, Descriptor> *ConvectiveNoForceRescaler<T, Descriptor>::clone() const
{
    return new ConvectiveNoForceRescaler<T, Descriptor>(*this);
}

/* ====================== Convective no force complete rescaler class ===================*/
template <typename T, template <typename U> class Descriptor>
ConvectiveNoForceCompleteRescaler<T, Descriptor>::ConvectiveNoForceCompleteRescaler() :
    Rescaler<T, Descriptor>()
{ }

template <typename T, template <typename U> class Descriptor>
void ConvectiveNoForceCompleteRescaler<T, Descriptor>::rescale(
    const Dynamics<T, Descriptor> &dyn, T xDt, std::vector<T> &rawData) const
{
    Array<T, Descriptor<T>::q> relFreq = dyn.getRelaxationFrequencies();

    if (Descriptor<T>::q + 1 + Descriptor<T>::d
        == ((rawData.size() - Descriptor<T>::ExternalField::numScalars)))
    {
        Array<T, Descriptor<T>::q> fNeq, mNeq;
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            fNeq[iPop] = rawData[1 + Descriptor<T>::d + iPop];
        }

        dynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::complete_ma2_moments(
            fNeq, mNeq);
        // rescale value of the relaxation time
        Array<T, Descriptor<T>::q> resRelFreq = this->computeRescaledRelFreq(relFreq, xDt);
        // only the non equilibrium moments are rescaled
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            T prefactor = relFreq[iPop] / resRelFreq[iPop] * xDt;
            // T prefactor = ((T)1-(T)0.5*relFreq[iPop]) / ((T)1-(T)0.5*resRelFreq[iPop]); // TODO
            // check: maybe in the wrong order
            mNeq[iPop] *= prefactor;
        }

        dynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::complete_ma2_populations(
            mNeq, fNeq);

        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            rawData[1 + Descriptor<T>::d + iPop] = fNeq[iPop];
        }
    } else {
        Array<T, Descriptor<T>::q> resRelFreq = this->computeRescaledRelFreq(relFreq, xDt);
        for (plint iPop = 0; iPop < SymmetricTensorImpl<T, Descriptor<T>::d>::n; ++iPop) {
            plint iA = 1 + Descriptor<T>::d + iPop;
            T prefactor = relFreq[iA] / resRelFreq[iA] * xDt;
            rawData[iA] *= prefactor;
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ConvectiveNoForceCompleteRescaler<T, Descriptor>
    *ConvectiveNoForceCompleteRescaler<T, Descriptor>::clone() const
{
    return new ConvectiveNoForceCompleteRescaler<T, Descriptor>(*this);
}

/* ====================== Convective no force complete rescaler class ===================*/
template <typename T, template <typename U> class Descriptor>
ConvectiveNoForceCompleteSmagorinskyRescaler<
    T, Descriptor>::ConvectiveNoForceCompleteSmagorinskyRescaler() :
    Rescaler<T, Descriptor>()
{ }

template <typename T, template <typename U> class Descriptor>
void ConvectiveNoForceCompleteSmagorinskyRescaler<T, Descriptor>::rescale(
    const Dynamics<T, Descriptor> &dyn, T xDt, std::vector<T> &rawData) const
{
    // PLB_ASSERT(false && "bug discovered in the convective no force smagorinsky rescaler. being
    // corrected.") PLB_ASSERT(order == 0 && "not implemented correctly for order == 1.")
    Array<T, Descriptor<T>::q> relFreq = dyn.getRelaxationFrequencies();
    Array<T, Descriptor<T>::q> fNeq, mNeq;
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        fNeq[iPop] = rawData[1 + Descriptor<T>::d + iPop];
    }

    dynamicsTemplatesImpl<T, typename Descriptor<T>::BaseDescriptor>::complete_ma2_moments(
        fNeq, mNeq);
    // rescale value of the relaxation time
    Array<T, Descriptor<T>::q> resRelFreq = this->computeRescaledRelFreq(relFreq, xDt);
    // only the non equilibrium moments are rescaled
    // for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
    //     mNeq[iPop] *= relFreqFactor;
    // }

    Array<T, SymmetricTensor<T, Descriptor>::n> piNeq;
    for (plint iPop = 1 + Descriptor<T>::d;
         iPop < 1 + Descriptor<T>::d + SymmetricTensor<T, Descriptor>::n; ++iPop)
    {
        piNeq[iPop - (1 + Descriptor<T>::d)] = mNeq[iPop];
    }

    T cSmago = dyn.getParameter(dynamicParams::smagorinskyConstant);
    Array<T, Descriptor<T>::q> sgsTerm = SmagoOperations<T, Descriptor>::computeSgsTensorTerm(
        Descriptor<T>::fullRho(rawData[0]), piNeq, cSmago,
        dyn.getParameter(dynamicParams::omega_shear));

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        T relFreqFactor = ((T)1 - (T)0.5 * relFreq[iPop]);  // TODO check: maybe in the wrong order
        T resRelFreqFactor = ((T)1 - (T)0.5 * resRelFreq[iPop]);
        T combFactor = (T)0.5 * (relFreq[iPop] - resRelFreq[iPop]);

        fNeq[iPop] *= relFreqFactor;
        fNeq[iPop] += combFactor * sgsTerm[iPop];
        fNeq[iPop] /= resRelFreqFactor;
    }

    // dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor >::complete_ma2_populations(
    // mNeq, fNeq );

    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        rawData[1 + Descriptor<T>::d + iPop] = fNeq[iPop];
    }
}

template <typename T, template <typename U> class Descriptor>
ConvectiveNoForceCompleteSmagorinskyRescaler<T, Descriptor>
    *ConvectiveNoForceCompleteSmagorinskyRescaler<T, Descriptor>::clone() const
{
    return new ConvectiveNoForceCompleteSmagorinskyRescaler<T, Descriptor>(*this);
}

}  // namespace plb

#endif  // GRID_REFINEMENT_HELPERS_HH
