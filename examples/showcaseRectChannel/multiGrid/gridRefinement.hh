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
#ifndef GRID_REFINEMENT_HH
#define GRID_REFINEMENT_HH

#include "multiGrid/gridRefinement.h"

namespace plb {

/** This function interpolates a function at the point (xValue,yValue)
 *  using the 2D Lagrangien formula using the vectors knownX, knownY and
 *  knownF. knownF must have knownX.size()*knownY.size() elements.
 */
template <typename T>
T interpolateValue(
    std::vector<T> const &knownX, std::vector<T> const &knownY,
    std::vector<std::vector<T> > const &knownF, T xValue, T yValue)
{
    // computation of the x coefficients of the polynomial
    std::vector<T> Linx;
    for (plint iX = 0; iX < (plint)knownX.size(); ++iX) {
        T res = 1.;
        for (plint i = 0; i < (plint)knownX.size(); ++i) {
            if (i != iX) {
                res *= (xValue - knownX[i]) / (knownX[iX] - knownX[i]);
            }
        }
        Linx.push_back(res);
    }

    // computation of the y coefficients of the polynomial
    std::vector<T> Limy;
    for (plint iY = 0; iY < (plint)knownY.size(); ++iY) {
        T res = 1.;
        for (plint i = 0; i < (plint)knownY.size(); ++i) {
            if (i != iY) {
                res *= (yValue - knownY[i]) / (knownY[iY] - knownY[i]);
            }
        }
        Limy.push_back(res);
    }

    // joining the whole following the given formula
    T result = 0.;
    for (plint iX = 0; iX < (plint)knownX.size(); ++iX) {
        for (plint iY = 0; iY < (plint)knownY.size(); ++iY) {
            result += Linx[iX] * Limy[iY] * knownF[iX][iY];
        }
    }
    return result;
}

/* *************** Class ConvectiveRescaleEngine **************************** */

template <typename T, template <typename U> class Descriptor>
const T ConvectiveRescaleEngine<T, Descriptor>::toFine_xDxInv = (T)2.;

template <typename T, template <typename U> class Descriptor>
const T ConvectiveRescaleEngine<T, Descriptor>::toFine_xDt = (T)1. / (T)2.;

template <typename T, template <typename U> class Descriptor>
const T ConvectiveRescaleEngine<T, Descriptor>::toCoarse_xDxInv = (T)1. / (T)2.;

template <typename T, template <typename U> class Descriptor>
const T ConvectiveRescaleEngine<T, Descriptor>::toCoarse_xDt = (T)2.;

template <typename T, template <typename U> class Descriptor>
ConvectiveRescaleEngine<T, Descriptor>::ConvectiveRescaleEngine(plint order_) : order(order_)
{ }

template <typename T, template <typename U> class Descriptor>
void ConvectiveRescaleEngine<T, Descriptor>::scaleCoarseFine(
    Cell<T, Descriptor> const &coarseCell, std::vector<T> &decomposedFineValues) const
{
    coarseCell.getDynamics().decompose(coarseCell, decomposedFineValues, order);

    if (!coarseCell.getDynamics().isAdvectionDiffusion()) {
        T tauCoarse = 1. / coarseCell.getDynamics().getOmega();
        T tauFine = 2. * tauCoarse - (T)1. / (T)2.;

        if (order == 0) {
            for (plint iFneq = 0; iFneq < Descriptor<T>::q; ++iFneq) {
                decomposedFineValues[1 + Descriptor<T>::d + iFneq] *= tauFine / tauCoarse;
            }
        } else {
            for (plint iPineq = 0; iPineq < SymmetricTensor<T, Descriptor>::n; ++iPineq) {
                decomposedFineValues[1 + Descriptor<T>::d + iPineq] *= tauFine / tauCoarse;
            }
        }

        coarseCell.getDynamics().rescale(decomposedFineValues, toFine_xDxInv, toFine_xDt, order);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConvectiveRescaleEngine<T, Descriptor>::scaleFineCoarse(
    Cell<T, Descriptor> const &fineCell, std::vector<T> &decomposedCoarseValues) const
{
    // decompose the fine cell populations
    fineCell.getDynamics().decompose(fineCell, decomposedCoarseValues, order);

    if (!fineCell.getDynamics().isAdvectionDiffusion()) {
        // computation of the factor associated with both tau
        T tauFine = 1. / fineCell.getDynamics().getOmega();
        T tauCoarse = (tauFine + (T)1. / (T)2.) / (T)2.;

        if (order == 0) {
            for (plint iFneq = 0; iFneq < Descriptor<T>::q; ++iFneq) {
                decomposedCoarseValues[1 + Descriptor<T>::d + iFneq] *= tauCoarse / tauFine;
            }
        } else {
            for (plint iPineq = 0; iPineq < SymmetricTensor<T, Descriptor>::n; ++iPineq) {
                decomposedCoarseValues[1 + Descriptor<T>::d + iPineq] *= tauCoarse / tauFine;
            }
        }

        // rescale and copy the values inside decomposedCoarseValues
        fineCell.getDynamics().rescale(
            decomposedCoarseValues, toCoarse_xDxInv, toCoarse_xDt, order);
    }
}

template <typename T, template <typename U> class Descriptor>
void ConvectiveRescaleEngine<T, Descriptor>::recompose(
    Cell<T, Descriptor> &cell, std::vector<T> const &decomposedValues) const
{
    cell.getDynamics().recompose(cell, decomposedValues, order);
}

template <typename T, template <typename U> class Descriptor>
plint ConvectiveRescaleEngine<T, Descriptor>::getDecompositionOrder() const
{
    return order;
}

template <typename T, template <typename U> class Descriptor>
ConvectiveRescaleEngine<T, Descriptor> *ConvectiveRescaleEngine<T, Descriptor>::clone() const
{
    return new ConvectiveRescaleEngine<T, Descriptor>(*this);
}

/* *************** Class NoScalingEngine ************************************ */

template <typename T, template <typename U> class Descriptor>
NoScalingEngine<T, Descriptor>::NoScalingEngine(plint order_) : order(order_)
{ }

template <typename T, template <typename U> class Descriptor>
void NoScalingEngine<T, Descriptor>::scaleCoarseFine(
    Cell<T, Descriptor> const &coarseCell, std::vector<T> &decomposedFineValues) const
{
    coarseCell.getDynamics().decompose(coarseCell, decomposedFineValues, order);
}

template <typename T, template <typename U> class Descriptor>
void NoScalingEngine<T, Descriptor>::scaleFineCoarse(
    Cell<T, Descriptor> const &fineCell, std::vector<T> &decomposedCoarseValues) const
{
    fineCell.getDynamics().decompose(fineCell, decomposedCoarseValues, order);
}

template <typename T, template <typename U> class Descriptor>
void NoScalingEngine<T, Descriptor>::recompose(
    Cell<T, Descriptor> &cell, std::vector<T> const &decomposedValues) const
{
    cell.getDynamics().recompose(cell, decomposedValues, order);
}

template <typename T, template <typename U> class Descriptor>
plint NoScalingEngine<T, Descriptor>::getDecompositionOrder() const
{
    return order;
}

template <typename T, template <typename U> class Descriptor>
NoScalingEngine<T, Descriptor> *NoScalingEngine<T, Descriptor>::clone() const
{
    return new NoScalingEngine<T, Descriptor>(*this);
}

/* ************* Interpolations ****************** */

/// [pop1  x  pop2] we interpolate over the x
template <typename T>
void linearInterpolation(
    std::vector<T> &pop1, std::vector<T> &pop2, std::vector<T> &decomposedValues)
{
    PLB_PRECONDITION(pop1.size() == pop2.size());
    decomposedValues.resize(pop1.size());

    for (pluint iVal = 0; iVal < decomposedValues.size(); ++iVal) {
        decomposedValues[iVal] = 1. / 2. * (pop1[iVal] + pop2[iVal]);
    }
}

/// [pop1 pop2 x pop3  pop4] we interpolate over the x
template <typename T>
void cubicCenteredInterpolation(
    std::vector<T> &pop1, std::vector<T> &pop2, std::vector<T> &pop3, std::vector<T> &pop4,
    std::vector<T> &decomposedValues)
{
    PLB_PRECONDITION(pop1.size() == pop2.size());
    decomposedValues.resize(pop1.size());

    for (pluint iVal = 0; iVal < decomposedValues.size(); ++iVal) {
        decomposedValues[iVal] =
            9. / 16. * (pop2[iVal] + pop3[iVal]) - 1. / 16. * (pop1[iVal] + pop4[iVal]);
    }
}

/// [pop1 x pop2 pop3] we interpolate over the x
template <typename T>
void quadraticNonCenteredInterpolation(
    std::vector<T> &pop1, std::vector<T> &pop2, std::vector<T> &pop3,
    std::vector<T> &decomposedValues)
{
    PLB_PRECONDITION(pop1.size() == pop2.size() && pop3.size() == pop1.size());
    decomposedValues.resize(pop1.size());

    for (pluint iVal = 0; iVal < decomposedValues.size(); ++iVal) {
        decomposedValues[iVal] = 3. / 8. * pop1[iVal] + 3. / 4. * pop2[iVal] - 1. / 8. * pop3[iVal];
    }
}

}  // namespace plb

#endif  // GRID_REFINEMENT_HH
