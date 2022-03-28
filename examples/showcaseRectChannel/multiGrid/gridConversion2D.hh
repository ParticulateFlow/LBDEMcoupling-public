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

#ifndef GRID_CONVERSION_2D_HH
#define GRID_CONVERSION_2D_HH

#include <algorithm>

#include "multiGrid/gridConversion2D.h"

namespace plb {

/* ******************* ScalarCopyFineToCoarseZerothOrder2D ******************* */

template <typename T>
ScalarCopyFineToCoarseZerothOrder2D<T>::ScalarCopyFineToCoarseZerothOrder2D(
    plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_) :
    dimDx(dimDx_), dimDt(dimDt_), dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be positive because the conversion goes from fine to coarse.
    PLB_ASSERT(dxScale >= 1);
}

template <typename T>
ScalarCopyFineToCoarseZerothOrder2D<T>::ScalarCopyFineToCoarseZerothOrder2D(
    ScalarCopyFineToCoarseZerothOrder2D<T> const &rhs) :
    dimDx(rhs.dimDx), dimDt(rhs.dimDt), dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T>
ScalarCopyFineToCoarseZerothOrder2D<T> &ScalarCopyFineToCoarseZerothOrder2D<T>::operator=(
    ScalarCopyFineToCoarseZerothOrder2D<T> const &rhs)
{
    dimDx = rhs.dimDx;
    dimDt = rhs.dimDt;
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T>
void ScalarCopyFineToCoarseZerothOrder2D<T>::process(
    Box2D fineDomain, ScalarField2D<T> &fineField, ScalarField2D<T> &coarseField)
{
    Dot2D posFine = fineField.getLocation();
    Dot2D posCoarse = coarseField.getLocation();

    plint stretch = util::twoToThePower(dxScale);
    Box2D coarseDomain(
        fineDomain.shift(posFine.x, posFine.y)
            .  // Convert to absolute fine coordinates.
        divideAndFitSmaller(stretch)
            .                                // Rescale, but don't exceed original domain.
        shift(-posCoarse.x, -posCoarse.y));  // Convert to relative coarse coordinates.
    PLB_ASSERT(contained(coarseDomain, coarseField.getBoundingBox()));

    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);

    plint fineX = (coarseDomain.x0 + posCoarse.x) * stretch - posFine.x;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX, fineX += stretch) {
        plint fineY = (coarseDomain.y0 + posCoarse.y) * stretch - posFine.y;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1;
             ++coarseY, fineY += stretch) {
            coarseField.get(coarseX, coarseY) = fineField.get(fineX, fineY) * scaleFactor;
        }
    }
}

template <typename T>
ScalarCopyFineToCoarseZerothOrder2D<T> *ScalarCopyFineToCoarseZerothOrder2D<T>::clone() const
{
    return new ScalarCopyFineToCoarseZerothOrder2D<T>(*this);
}

/* ******************* TensorCopyFineToCoarseZerothOrder2D ******************* */

template <typename T, int nDim>
TensorCopyFineToCoarseZerothOrder2D<T, nDim>::TensorCopyFineToCoarseZerothOrder2D(
    plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_) :
    dimDx(dimDx_), dimDt(dimDt_), dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be positive because the conversion goes from fine to coarse.
    PLB_ASSERT(dxScale >= 1);
}

template <typename T, int nDim>
TensorCopyFineToCoarseZerothOrder2D<T, nDim>::TensorCopyFineToCoarseZerothOrder2D(
    TensorCopyFineToCoarseZerothOrder2D<T, nDim> const &rhs) :
    dimDx(rhs.dimDx), dimDt(rhs.dimDt), dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T, int nDim>
TensorCopyFineToCoarseZerothOrder2D<T, nDim>
    &TensorCopyFineToCoarseZerothOrder2D<T, nDim>::operator=(
        TensorCopyFineToCoarseZerothOrder2D<T, nDim> const &rhs)
{
    dimDx = rhs.dimDx;
    dimDt = rhs.dimDt;
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T, int nDim>
void TensorCopyFineToCoarseZerothOrder2D<T, nDim>::process(
    Box2D fineDomain, TensorField2D<T, nDim> &fineField, TensorField2D<T, nDim> &coarseField)
{
    Dot2D posFine = fineField.getLocation();
    Dot2D posCoarse = coarseField.getLocation();

    plint stretch = util::twoToThePower(dxScale);
    Box2D coarseDomain(
        fineDomain.shift(posFine.x, posFine.y)
            .  // Convert to absolute fine coordinates.
        divideAndFitSmaller(stretch)
            .                                // Rescale, but don't exceed original domain.
        shift(-posCoarse.x, -posCoarse.y));  // Convert to relative coarse coordinates.
    PLB_ASSERT(contained(coarseDomain, coarseField.getBoundingBox()));

    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);

    plint fineX = (coarseDomain.x0 + posCoarse.x) * stretch - posFine.x;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX, fineX += stretch) {
        plint fineY = (coarseDomain.y0 + posCoarse.y) * stretch - posFine.y;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1;
             ++coarseY, fineY += stretch) {
            coarseField.get(coarseX, coarseY) = fineField.get(fineX, fineY) * scaleFactor;
        }
    }
}

template <typename T, int nDim>
TensorCopyFineToCoarseZerothOrder2D<T, nDim> *TensorCopyFineToCoarseZerothOrder2D<T, nDim>::clone()
    const
{
    return new TensorCopyFineToCoarseZerothOrder2D<T, nDim>(*this);
}

/* ******************* LatticeCopyFineToCoarseZerothOrder2D ******************* */

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor>::LatticeCopyFineToCoarseZerothOrder2D(
    plint dxScale_, plint dtScale_) :
    dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be positive because the conversion goes from fine to coarse.
    PLB_ASSERT(dxScale >= 1);
}

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor>::LatticeCopyFineToCoarseZerothOrder2D(
    LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor> const &rhs) :
    dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor>
    &LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor>::operator=(
        LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor> const &rhs)
{
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor>::process(
    Box2D fineDomain, BlockLattice2D<T, Descriptor> &fineLattice,
    BlockLattice2D<T, Descriptor> &coarseLattice)
{
    static const int velDimDx = 1;
    static const int velDimDt = -1;
    static const int presDimDx = 2;
    static const int presDimDt = -2;
    static const int viscDimDx = 2;
    static const int viscDimDt = -1;
    static const int fNeqDimDx = 0;
    static const int fNeqDimDt = -1;
    T velScale = scaleFromReference(dxScale, velDimDx, dtScale, velDimDt);
    T presScale = scaleFromReference(dxScale, presDimDx, dtScale, presDimDt);
    T viscScale = scaleFromReference(dxScale, viscDimDx, dtScale, viscDimDt);
    T fNeqScale = scaleFromReference(dxScale, fNeqDimDx, dtScale, fNeqDimDt);

    Dot2D posFine = fineLattice.getLocation();
    Dot2D posCoarse = coarseLattice.getLocation();

    plint stretch = util::twoToThePower(dxScale);
    Box2D coarseDomain(
        fineDomain.shift(posFine.x, posFine.y)
            .  // Convert to absolute fine coordinates.
        divideAndFitSmaller(stretch)
            .                                // Rescale, but don't exceed original domain.
        shift(-posCoarse.x, -posCoarse.y));  // Convert to relative coarse coordinates.
    PLB_ASSERT(contained(coarseDomain, coarseLattice.getBoundingBox()));

    plint fineX = (coarseDomain.x0 + posCoarse.x) * stretch - posFine.x;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX, fineX += stretch) {
        plint fineY = (coarseDomain.y0 + posCoarse.y) * stretch - posFine.y;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1;
             ++coarseY, fineY += stretch) {
            // Fine --> Coarse
            Cell<T, Descriptor> const &fineCell = fineLattice.get(fineX, fineY);
            coarseLattice.attributeDynamics(coarseX, coarseY, fineCell.getDynamics().clone());
            Cell<T, Descriptor> &coarseCell = coarseLattice.get(coarseX, coarseY);
            T fineOmega = fineCell.getDynamics().getOmega();
            T fine_nu_cs2 = 1. / fineOmega - 0.5;
            T coarse_nu_cs2 = viscScale * fine_nu_cs2;
            T coarseOmega = 1. / (coarse_nu_cs2 + 0.5);
            coarseCell.getDynamics().setOmega(coarseOmega);

            T localFneqScale = fNeqScale * fineOmega / coarseOmega;
            std::vector<T> components;
            plint order = 0;
            fineCell.getDynamics().decompose(fineCell, components, order);
            T rhoFine = Descriptor<T>::fullRho(components[0]);
            if (!fineCell.getDynamics().velIsJ()) {
                for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
                    components[1 + iDim] /= rhoFine;
                }
            }
            if (components.size() == 1 + Descriptor<T>::d + Descriptor<T>::q) {
                // PLB_ASSERT(components.size()==1+Descriptor<T>::d+Descriptor<T>::q);
                // 1. Density.
                T p_cs2 = components[0] - (1 - Descriptor<T>::SkordosFactor());
                components[0] = (1 - Descriptor<T>::SkordosFactor()) + p_cs2 * presScale;
                // 2. Velocity.
                for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
                    components[1 + iDim] *= velScale;
                }
                // 3. Off-equilibrium populations.
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    components[1 + Descriptor<T>::d + iPop] *= localFneqScale;
                }
                T rhoCoarse = Descriptor<T>::fullRho(components[0]);
                if (!coarseCell.getDynamics().velIsJ()) {
                    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
                        components[1 + iDim] *= rhoCoarse;
                    }
                }
            }
            coarseCell.getDynamics().recompose(coarseCell, components, order);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor>
    *LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor>::clone() const
{
    return new LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor>(*this);
}

/* ******************* ScalarInterpolateCoarseToFine2D ******************* */

template <typename T>
ScalarInterpolateCoarseToFine2D<T>::ScalarInterpolateCoarseToFine2D(
    plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_) :
    dimDx(dimDx_), dimDt(dimDt_), dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be negative because the conversion goes from coarse to fine.
    PLB_ASSERT(dxScale <= -1);
}

template <typename T>
ScalarInterpolateCoarseToFine2D<T>::ScalarInterpolateCoarseToFine2D(
    ScalarInterpolateCoarseToFine2D<T> const &rhs) :
    dimDx(rhs.dimDx), dimDt(rhs.dimDt), dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T>
ScalarInterpolateCoarseToFine2D<T> &ScalarInterpolateCoarseToFine2D<T>::operator=(
    ScalarInterpolateCoarseToFine2D<T> const &rhs)
{
    dimDx = rhs.dimDx;
    dimDt = rhs.dimDt;
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T>
void ScalarInterpolateCoarseToFine2D<T>::process(
    Box2D coarseDomain, ScalarField2D<T> &coarseField, ScalarField2D<T> &fineField)
{
    Dot2D posFine = fineField.getLocation();
    Dot2D posCoarse = coarseField.getLocation();

    plint stretch = util::twoToThePower(-dxScale);
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);

    std::vector<Dot2D> cellPos;
    std::vector<T> weights;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX) {
        plint fineX = (coarseX + posCoarse.x) * stretch - posFine.x;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1; ++coarseY) {
            plint fineY = (coarseY + posCoarse.y) * stretch - posFine.y;
            for (plint iX = fineX - stretch / 2 + 1; iX <= fineX + stretch / 2; ++iX) {
                for (plint iY = fineY - stretch / 2 + 1; iY <= fineY + stretch / 2; ++iY) {
                    if (contained(iX, iY, fineField.getBoundingBox())) {
                        T coarseIx = (T)(iX + posFine.x) / (T)stretch;
                        T coarseIy = (T)(iY + posFine.y) / (T)stretch;
                        linearInterpolationCoefficients(
                            coarseField, Array<T, 2>(coarseIx, coarseIy), cellPos, weights);
                        fineField.get(iX, iY) = 0.;
                        for (plint i = 0; i < 4; ++i) {
                            fineField.get(iX, iY) +=
                                weights[i] * coarseField.get(cellPos[i].x, cellPos[i].y);
                        }
                        fineField.get(iX, iY) *= scaleFactor;
                    }
                }
            }
        }
    }
}

template <typename T>
ScalarInterpolateCoarseToFine2D<T> *ScalarInterpolateCoarseToFine2D<T>::clone() const
{
    return new ScalarInterpolateCoarseToFine2D<T>(*this);
}

/* ******************* TensorInterpolateCoarseToFine2D ******************* */

template <typename T, int nDim>
TensorInterpolateCoarseToFine2D<T, nDim>::TensorInterpolateCoarseToFine2D(
    plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_) :
    dimDx(dimDx_), dimDt(dimDt_), dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be negative because the conversion goes from coarse to fine.
    PLB_ASSERT(dxScale <= -1);
}

template <typename T, int nDim>
TensorInterpolateCoarseToFine2D<T, nDim>::TensorInterpolateCoarseToFine2D(
    TensorInterpolateCoarseToFine2D<T, nDim> const &rhs) :
    dimDx(rhs.dimDx), dimDt(rhs.dimDt), dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T, int nDim>
TensorInterpolateCoarseToFine2D<T, nDim> &TensorInterpolateCoarseToFine2D<T, nDim>::operator=(
    TensorInterpolateCoarseToFine2D<T, nDim> const &rhs)
{
    dimDx = rhs.dimDx;
    dimDt = rhs.dimDt;
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T, int nDim>
void TensorInterpolateCoarseToFine2D<T, nDim>::process(
    Box2D coarseDomain, TensorField2D<T, nDim> &coarseField, TensorField2D<T, nDim> &fineField)
{
    Dot2D posFine = fineField.getLocation();
    Dot2D posCoarse = coarseField.getLocation();

    plint stretch = util::twoToThePower(-dxScale);
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);

    std::vector<Dot2D> cellPos;
    std::vector<T> weights;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX) {
        plint fineX = (coarseX + posCoarse.x) * stretch - posFine.x;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1; ++coarseY) {
            plint fineY = (coarseY + posCoarse.y) * stretch - posFine.y;
            for (plint iX = fineX - stretch / 2 + 1; iX <= fineX + stretch / 2; ++iX) {
                for (plint iY = fineY - stretch / 2 + 1; iY <= fineY + stretch / 2; ++iY) {
                    if (contained(iX, iY, fineField.getBoundingBox())) {
                        T coarseIx = (T)(iX + posFine.x) / (T)stretch;
                        T coarseIy = (T)(iY + posFine.y) / (T)stretch;
                        linearInterpolationCoefficients(
                            coarseField, Array<T, 2>(coarseIx, coarseIy), cellPos, weights);
                        fineField.get(iX, iY).resetToZero();
                        for (plint i = 0; i < 4; ++i) {
                            fineField.get(iX, iY) +=
                                weights[i] * coarseField.get(cellPos[i].x, cellPos[i].y);
                        }
                        fineField.get(iX, iY) *= scaleFactor;
                    }
                }
            }
        }
    }
}

template <typename T, int nDim>
TensorInterpolateCoarseToFine2D<T, nDim> *TensorInterpolateCoarseToFine2D<T, nDim>::clone() const
{
    return new TensorInterpolateCoarseToFine2D<T, nDim>(*this);
}

/* ******************* LatticeInterpolateCoarseToFine2D ******************* */

template <typename T, template <typename U> class Descriptor>
LatticeInterpolateCoarseToFine2D<T, Descriptor>::LatticeInterpolateCoarseToFine2D(
    plint dxScale_, plint dtScale_) :
    dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be negative because the conversion goes from coarse to fine.
    PLB_ASSERT(dxScale <= -1);
}

template <typename T, template <typename U> class Descriptor>
LatticeInterpolateCoarseToFine2D<T, Descriptor>::LatticeInterpolateCoarseToFine2D(
    LatticeInterpolateCoarseToFine2D<T, Descriptor> const &rhs) :
    dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T, template <typename U> class Descriptor>
LatticeInterpolateCoarseToFine2D<T, Descriptor>
    &LatticeInterpolateCoarseToFine2D<T, Descriptor>::operator=(
        LatticeInterpolateCoarseToFine2D<T, Descriptor> const &rhs)
{
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void LatticeInterpolateCoarseToFine2D<T, Descriptor>::process(
    Box2D coarseDomain, BlockLattice2D<T, Descriptor> &coarseLattice,
    BlockLattice2D<T, Descriptor> &fineLattice)
{
    static const int velDimDx = 1;
    static const int velDimDt = -1;
    static const int presDimDx = 2;
    static const int presDimDt = -2;
    // static const int viscDimDx =  2;
    // static const int viscDimDt = -1;
    static const int fNeqDimDx = 0;
    static const int fNeqDimDt = -1;
    T velScale = scaleFromReference(dxScale, velDimDx, dtScale, velDimDt);
    T presScale = scaleFromReference(dxScale, presDimDx, dtScale, presDimDt);
    // T viscScale = scaleFromReference(dxScale, viscDimDx, dtScale, viscDimDt);
    T fNeqScale = scaleFromReference(dxScale, fNeqDimDx, dtScale, fNeqDimDt);

    Dot2D posFine = fineLattice.getLocation();
    Dot2D posCoarse = coarseLattice.getLocation();

    plint stretch = util::intTwoToThePower(-dxScale);

    std::vector<Dot2D> cellPos;
    std::vector<T> weights;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX) {
        plint fineX = (coarseX + posCoarse.x) * stretch - posFine.x;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1; ++coarseY) {
            plint fineY = (coarseY + posCoarse.y) * stretch - posFine.y;
            for (plint iX = fineX - stretch / 2 + 1; iX <= fineX + stretch / 2; ++iX) {
                for (plint iY = fineY - stretch / 2 + 1; iY <= fineY + stretch / 2; ++iY) {
                    if (contained(iX, iY, fineLattice.getBoundingBox())) {
                        T coarseIx = (T)(iX + posFine.x) / (T)stretch;
                        T coarseIy = (T)(iY + posFine.y) / (T)stretch;
                        linearInterpolationCoefficients(
                            coarseLattice, Array<T, 2>(coarseIx, coarseIy), cellPos, weights);
                        plint order = 0;
                        plint numComponents = 1 + Descriptor<T>::d + Descriptor<T>::q
                                              + Descriptor<T>::ExternalField::numScalars;
                        // TODO: For now external scalars are set to zero during conversion. They
                        // should be
                        //   rescaled properly.
                        std::vector<T> components(numComponents), oldComponents(numComponents);
                        std::fill(components.begin(), components.end(), 0.);
                        // T coarse_nu_cs2 = 0.;
                        components[0] = 1 - Descriptor<T>::SkordosFactor();  // Density offset.
                        for (plint i = 0; i < 4; ++i) {
                            Cell<T, Descriptor> const &coarseCell =
                                coarseLattice.get(cellPos[i].x, cellPos[i].y);
                            // T coarseOmega = coarseCell.getDynamics().getOmega();
                            // T tmp_coarse_nu_cs2 = 1./coarseOmega - 0.5;
                            // coarse_nu_cs2 += weights[i]*tmp_coarse_nu_cs2;
                            coarseCell.getDynamics().decompose(coarseCell, oldComponents, order);
                            // 1. rhoBar. Keep the constant term (1), and rescale the remaining term
                            // like a pressure.
                            T pres_cs2 = oldComponents[0] - (1 - Descriptor<T>::SkordosFactor());
                            components[0] += weights[i] * pres_cs2 * presScale;
                            // 2. j.
                            T velMultFactor = weights[i] * velScale;
                            for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
                                components[1 + iDim] += velMultFactor * oldComponents[1 + iDim];
                            }
                            // 3. Off-equilibrium populations.
                            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                                components[1 + Descriptor<T>::d + iPop] +=
                                    weights[i] * fNeqScale
                                    * oldComponents[1 + Descriptor<T>::d + iPop];
                            }
                        }
                        Cell<T, Descriptor> &fineCell = fineLattice.get(iX, iY);
                        // Rescaling of the relaxation time has been deactivated for now, because
                        // it produces numerical instability when a checkpoint file is read and
                        // scaled. Needs to be fixed at some point.
                        //
                        // T fine_nu_cs2 = coarse_nu_cs2*viscScale;
                        // T coarseOmega = 1. / (coarse_nu_cs2+0.5);
                        // T fineOmega = 1. / (fine_nu_cs2+0.5);
                        // T omegaRatio = coarseOmega / fineOmega;
                        // for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                        //    components[1+Descriptor<T>::d+iPop] *= omegaRatio;
                        //}
                        // fineCell.getDynamics().setOmega(fineOmega);

                        fineCell.getDynamics().recompose(fineCell, components, order);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeInterpolateCoarseToFine2D<T, Descriptor>
    *LatticeInterpolateCoarseToFine2D<T, Descriptor>::clone() const
{
    return new LatticeInterpolateCoarseToFine2D<T, Descriptor>(*this);
}

/* ******************* LatticeLinearInterpolateCoarseToFine2D ******************* */

template <typename T, template <typename U> class Descriptor>
LatticeLinearInterpolateCoarseToFine2D<T, Descriptor>::LatticeLinearInterpolateCoarseToFine2D(
    plint dxScale_, plint dtScale_) :
    dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be negative because the conversion goes from coarse to fine.
    PLB_ASSERT(dxScale <= -1);
}

template <typename T, template <typename U> class Descriptor>
LatticeLinearInterpolateCoarseToFine2D<T, Descriptor>::LatticeLinearInterpolateCoarseToFine2D(
    LatticeLinearInterpolateCoarseToFine2D<T, Descriptor> const &rhs) :
    dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T, template <typename U> class Descriptor>
LatticeLinearInterpolateCoarseToFine2D<T, Descriptor>
    &LatticeLinearInterpolateCoarseToFine2D<T, Descriptor>::operator=(
        LatticeLinearInterpolateCoarseToFine2D<T, Descriptor> const &rhs)
{
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void LatticeLinearInterpolateCoarseToFine2D<T, Descriptor>::process(
    Box2D coarseDomain, BlockLattice2D<T, Descriptor> &coarseLattice,
    BlockLattice2D<T, Descriptor> &fineLattice)
{
    static const int velDimDx = 1;
    static const int velDimDt = -1;
    static const int presDimDx = 2;
    static const int presDimDt = -2;
    static const int viscDimDx = 2;
    static const int viscDimDt = -1;
    static const int fNeqDimDx = 0;
    static const int fNeqDimDt = -1;
    T velScale = scaleFromReference(dxScale, velDimDx, dtScale, velDimDt);
    T presScale = scaleFromReference(dxScale, presDimDx, dtScale, presDimDt);
    T viscScale = scaleFromReference(dxScale, viscDimDx, dtScale, viscDimDt);
    T fNeqScale = scaleFromReference(dxScale, fNeqDimDx, dtScale, fNeqDimDt);

    Dot2D posFine = fineLattice.getLocation();
    Dot2D posCoarse = coarseLattice.getLocation();

    plint stretch = util::intTwoToThePower(-dxScale);

    std::vector<Dot2D> cellPos;
    std::vector<T> weights;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX) {
        plint fineX = (coarseX + posCoarse.x) * stretch - posFine.x;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1; ++coarseY) {
            plint fineY = (coarseY + posCoarse.y) * stretch - posFine.y;
            for (plint iX = fineX - stretch / 2 + 1; iX <= fineX + stretch / 2; ++iX) {
                for (plint iY = fineY - stretch / 2 + 1; iY <= fineY + stretch / 2; ++iY) {
                    if (contained(iX, iY, fineLattice.getBoundingBox())) {
                        T coarseIx = (T)(iX + posFine.x) / (T)stretch;
                        T coarseIy = (T)(iY + posFine.y) / (T)stretch;
                        linearInterpolationCoefficients(
                            coarseLattice, Array<T, 2>(coarseIx, coarseIy), cellPos, weights);
                        plint order = 0;
                        plint numComponents = 1 + Descriptor<T>::d + Descriptor<T>::q
                                              + Descriptor<T>::ExternalField::numScalars;
                        // TODO: For now external scalars are set to zero during conversion. They
                        // should be
                        //   rescaled properly.
                        std::vector<T> components(numComponents), oldComponents(numComponents);
                        std::fill(components.begin(), components.end(), 0.);
                        T coarse_nu_cs2 = 0.;
                        components[0] = 1 - Descriptor<T>::SkordosFactor();  // Density offset.
                        for (plint i = 0; i < 4; ++i) {
                            Cell<T, Descriptor> const &coarseCell =
                                coarseLattice.get(cellPos[i].x, cellPos[i].y);
                            T coarseOmega = coarseCell.getDynamics().getOmega();
                            T tmp_coarse_nu_cs2 = 1. / coarseOmega - 0.5;
                            coarse_nu_cs2 += weights[i] * tmp_coarse_nu_cs2;
                            coarseCell.getDynamics().decompose(coarseCell, oldComponents, order);
                            // 1. rhoBar. Keep the constant term (1), and rescale the remaining term
                            // like a pressure.
                            T pres_cs2 = oldComponents[0] - (1 - Descriptor<T>::SkordosFactor());
                            components[0] += weights[i] * pres_cs2 * presScale;
                            // 2. j.
                            T velMultFactor = weights[i] * velScale;
                            for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
                                components[1 + iDim] += velMultFactor * oldComponents[1 + iDim];
                            }
                            // 3. Off-equilibrium populations.
                            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                                components[1 + Descriptor<T>::d + iPop] +=
                                    weights[i] * fNeqScale
                                    * oldComponents[1 + Descriptor<T>::d + iPop];
                            }
                        }
                        Cell<T, Descriptor> &fineCell = fineLattice.get(iX, iY);

                        T fine_nu_cs2 = coarse_nu_cs2 * viscScale;
                        T coarseOmega = 1. / (coarse_nu_cs2 + 0.5);
                        T fineOmega = 1. / (fine_nu_cs2 + 0.5);
                        T omegaRatio = coarseOmega / fineOmega;
                        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                            components[1 + Descriptor<T>::d + iPop] *= omegaRatio;
                        }
                        fineCell.getDynamics().setOmega(fineOmega);

                        fineCell.getDynamics().recompose(fineCell, components, order);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeLinearInterpolateCoarseToFine2D<T, Descriptor>
    *LatticeLinearInterpolateCoarseToFine2D<T, Descriptor>::clone() const
{
    return new LatticeLinearInterpolateCoarseToFine2D<T, Descriptor>(*this);
}

/* ******************* Free functions ******************* */

inline Box2D scaleBox(Box2D box, plint nLevel)
{
    if (nLevel > 0) {
        return box.multiply(util::intTwoToThePower(nLevel));
    } else if (nLevel < 0) {
        // If the fine-grid box does not fit with the coarse-grid box,
        //   let's shrink the coarse-grid box by one fine cell. This
        //   makes sure we'll never exceed the fine-grid box
        return box.divideAndFitSmaller(util::intTwoToThePower(-nLevel));
    } else {
        return box;
    }
}

inline MultiBlockManagement2D scaleMultiBlockManagement(
    MultiBlockManagement2D const &multiBlockManagement, plint nLevel)
{
    SparseBlockStructure2D const &sparseBlock = multiBlockManagement.getSparseBlockStructure();

    Box2D rescaledBoundingBox = scaleBox(sparseBlock.getBoundingBox(), nLevel);
    SparseBlockStructure2D scaledSparseBlock(rescaledBoundingBox);
    plint stretch = 1;
    if (nLevel > 0) {
        stretch = util::intTwoToThePower(nLevel);
    }

    std::map<plint, Box2D>::const_iterator it = sparseBlock.getBulks().begin();
    for (; it != sparseBlock.getBulks().end(); ++it) {
        Box2D scaledBulk = scaleBox(it->second, nLevel);
        if (stretch > 1) {
            scaledBulk.x0 -= stretch / 2 - 1;
            scaledBulk.y0 -= stretch / 2 - 1;
            scaledBulk.x1 += stretch / 2;
            scaledBulk.y1 += stretch / 2;
        }
        Box2D scaledUniqueBulk = scaledBulk;  // TODO: compute unique bulk properly.
        plint blockId = it->first;
        scaledSparseBlock.addBlock(scaledBulk, scaledUniqueBulk, blockId);
    }
    return MultiBlockManagement2D(
        scaledSparseBlock, multiBlockManagement.getThreadAttribution().clone(),
        multiBlockManagement.getEnvelopeWidth(),
        multiBlockManagement.getRefinementLevel() + nLevel);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > coarsen(
    MultiScalarField2D<T> &fineField, plint dimDx, plint dimDt, plint dxScale, plint dtScale)
{
    PLB_PRECONDITION(dxScale >= 1);
    // Relative level is negative when going from fine to coarse.
    plint relativeLevel = -dxScale;
    MultiBlockManagement2D management =
        scaleMultiBlockManagement(fineField.getMultiBlockManagement(), relativeLevel);
    MultiScalarField2D<T> *result = new MultiScalarField2D<T>(
        management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>());
    applyProcessingFunctional(
        new ScalarCopyFineToCoarseZerothOrder2D<T>(dimDx, dimDt, dxScale, dtScale),
        fineField.getBoundingBox(), fineField, *result);

    return std::unique_ptr<MultiScalarField2D<T> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > coarsen(
    MultiTensorField2D<T, nDim> &fineField, plint dimDx, plint dimDt, plint dxScale, plint dtScale)
{
    PLB_PRECONDITION(dxScale >= 1);
    // Relative level is negative when going from fine to coarse.
    plint relativeLevel = -dxScale;
    MultiBlockManagement2D management =
        scaleMultiBlockManagement(fineField.getMultiBlockManagement(), relativeLevel);
    MultiTensorField2D<T, nDim> *result = new MultiTensorField2D<T, nDim>(
        management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>());
    applyProcessingFunctional(
        new TensorCopyFineToCoarseZerothOrder2D<T, nDim>(dimDx, dimDt, dxScale, dtScale),
        fineField.getBoundingBox(), fineField, *result);

    return std::unique_ptr<MultiTensorField2D<T, nDim> >(result);
}

template <typename T>
std::unique_ptr<MultiScalarField2D<T> > refine(
    MultiScalarField2D<T> &coarseField, plint dimDx, plint dimDt, plint dxScale, plint dtScale)
{
    PLB_PRECONDITION(dxScale <= -1);
    // Relative level is positive when going from coarse to fine.
    plint relativeLevel = -dxScale;
    MultiBlockManagement2D management =
        scaleMultiBlockManagement(coarseField.getMultiBlockManagement(), relativeLevel);
    MultiScalarField2D<T> *result = new MultiScalarField2D<T>(
        management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiScalarAccess<T>());
    applyProcessingFunctional(
        new ScalarInterpolateCoarseToFine2D<T>(dimDx, dimDt, dxScale, dtScale),
        coarseField.getBoundingBox(), coarseField, *result);

    return std::unique_ptr<MultiScalarField2D<T> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField2D<T, nDim> > refine(
    MultiTensorField2D<T, nDim> &coarseField, plint dimDx, plint dimDt, plint dxScale,
    plint dtScale)
{
    PLB_PRECONDITION(dxScale <= -1);
    // Relative level is positive when going from coarse to fine.
    plint relativeLevel = -dxScale;
    MultiBlockManagement2D management =
        scaleMultiBlockManagement(coarseField.getMultiBlockManagement(), relativeLevel);
    MultiTensorField2D<T, nDim> *result = new MultiTensorField2D<T, nDim>(
        management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>());
    applyProcessingFunctional(
        new TensorInterpolateCoarseToFine2D<T, nDim>(dimDx, dimDt, dxScale, dtScale),
        coarseField.getBoundingBox(), coarseField, *result);

    return std::unique_ptr<MultiTensorField2D<T, nDim> >(result);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > coarsen(
    MultiBlockLattice2D<T, Descriptor> &fineLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics)
{
    PLB_PRECONDITION(dxScale >= 1);
    // Relative level is negative when going from fine to coarse.
    plint relativeLevel = -dxScale;
    MultiBlockManagement2D management =
        scaleMultiBlockManagement(fineLattice.getMultiBlockManagement(), relativeLevel);
    MultiBlockLattice2D<T, Descriptor> *result = new MultiBlockLattice2D<T, Descriptor>(
        management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
        backgroundDynamics->clone());
    defineDynamics(*result, result->getBoundingBox(), backgroundDynamics);
    applyProcessingFunctional(
        new LatticeCopyFineToCoarseZerothOrder2D<T, Descriptor>(dxScale, dtScale),
        fineLattice.getBoundingBox(), fineLattice, *result);

    return std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(result);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > refine(
    MultiBlockLattice2D<T, Descriptor> &coarseLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics)
{
    PLB_PRECONDITION(dxScale <= -1);
    // Relative level is positive when going from coarse to fine.
    plint relativeLevel = -dxScale;
    MultiBlockManagement2D management =
        scaleMultiBlockManagement(coarseLattice.getMultiBlockManagement(), relativeLevel);
    MultiBlockLattice2D<T, Descriptor> *result = new MultiBlockLattice2D<T, Descriptor>(
        management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
        backgroundDynamics->clone());
    defineDynamics(*result, result->getBoundingBox(), backgroundDynamics);
    applyProcessingFunctional(
        new LatticeInterpolateCoarseToFine2D<T, Descriptor>(dxScale, dtScale),
        coarseLattice.getBoundingBox(), coarseLattice, *result);

    return std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(result);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > linearRefine(
    MultiBlockLattice2D<T, Descriptor> &coarseLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics)
{
    PLB_PRECONDITION(dxScale <= -1);
    // Relative level is positive when going from coarse to fine.
    plint relativeLevel = -dxScale;
    MultiBlockManagement2D management =
        scaleMultiBlockManagement(coarseLattice.getMultiBlockManagement(), relativeLevel);
    MultiBlockLattice2D<T, Descriptor> *result = new MultiBlockLattice2D<T, Descriptor>(
        management, defaultMultiBlockPolicy2D().getBlockCommunicator(),
        defaultMultiBlockPolicy2D().getCombinedStatistics(),
        defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
        backgroundDynamics->clone());
    defineDynamics(*result, result->getBoundingBox(), backgroundDynamics);
    applyProcessingFunctional(
        new LatticeLinearInterpolateCoarseToFine2D<T, Descriptor>(dxScale, dtScale),
        coarseLattice.getBoundingBox(), coarseLattice, *result);

    return std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(result);
}

}  // namespace plb

#endif  // GRID_CONVERSION_2D_HH
