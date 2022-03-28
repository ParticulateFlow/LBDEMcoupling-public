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

#ifndef GRID_CONVERSION_3D_HH
#define GRID_CONVERSION_3D_HH

#include <algorithm>

#include "finiteDifference/interpolations3D.h"
#include "multiGrid/gridConversion3D.h"

namespace plb {

/* ******************* ScalarCopyFineToCoarseZerothOrder3D ******************* */

template <typename T>
ScalarCopyFineToCoarseZerothOrder3D<T>::ScalarCopyFineToCoarseZerothOrder3D(
    plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_) :
    dimDx(dimDx_), dimDt(dimDt_), dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be positive because the conversion goes from fine to coarse.
    PLB_ASSERT(dxScale >= 1);
}

template <typename T>
ScalarCopyFineToCoarseZerothOrder3D<T>::ScalarCopyFineToCoarseZerothOrder3D(
    ScalarCopyFineToCoarseZerothOrder3D<T> const &rhs) :
    dimDx(rhs.dimDx), dimDt(rhs.dimDt), dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T>
ScalarCopyFineToCoarseZerothOrder3D<T> &ScalarCopyFineToCoarseZerothOrder3D<T>::operator=(
    ScalarCopyFineToCoarseZerothOrder3D<T> const &rhs)
{
    dimDx = rhs.dimDx;
    dimDt = rhs.dimDt;
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T>
void ScalarCopyFineToCoarseZerothOrder3D<T>::process(
    Box3D fineDomain, ScalarField3D<T> &fineField, ScalarField3D<T> &coarseField)
{
    Dot3D posFine = fineField.getLocation();
    Dot3D posCoarse = coarseField.getLocation();

    plint stretch = util::twoToThePower(dxScale);
    Box3D coarseDomain(
        fineDomain.shift(posFine.x, posFine.y, posFine.z)
            .  // Convert to absolute fine coordinates.
        divideAndFitSmaller(stretch)
            .  // Rescale, but don't exceed original domain.
        shift(
            -posCoarse.x, -posCoarse.y, -posCoarse.z));  // Convert to relative coarse coordinates.
    PLB_ASSERT(contained(coarseDomain, coarseField.getBoundingBox()));

    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);

    plint fineX = (coarseDomain.x0 + posCoarse.x) * stretch - posFine.x;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX, fineX += stretch) {
        plint fineY = (coarseDomain.y0 + posCoarse.y) * stretch - posFine.y;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1;
             ++coarseY, fineY += stretch) {
            plint fineZ = (coarseDomain.z0 + posCoarse.z) * stretch - posFine.z;
            for (plint coarseZ = coarseDomain.z0; coarseZ <= coarseDomain.z1;
                 ++coarseZ, fineZ += stretch) {
                coarseField.get(coarseX, coarseY, coarseZ) =
                    fineField.get(fineX, fineY, fineZ) * scaleFactor;
            }
        }
    }
}

template <typename T>
ScalarCopyFineToCoarseZerothOrder3D<T> *ScalarCopyFineToCoarseZerothOrder3D<T>::clone() const
{
    return new ScalarCopyFineToCoarseZerothOrder3D<T>(*this);
}

/* ******************* TensorCopyFineToCoarseZerothOrder3D ******************* */

template <typename T, int nDim>
TensorCopyFineToCoarseZerothOrder3D<T, nDim>::TensorCopyFineToCoarseZerothOrder3D(
    plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_) :
    dimDx(dimDx_), dimDt(dimDt_), dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be positive because the conversion goes from fine to coarse.
    PLB_ASSERT(dxScale >= 1);
}

template <typename T, int nDim>
TensorCopyFineToCoarseZerothOrder3D<T, nDim>::TensorCopyFineToCoarseZerothOrder3D(
    TensorCopyFineToCoarseZerothOrder3D<T, nDim> const &rhs) :
    dimDx(rhs.dimDx), dimDt(rhs.dimDt), dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T, int nDim>
TensorCopyFineToCoarseZerothOrder3D<T, nDim>
    &TensorCopyFineToCoarseZerothOrder3D<T, nDim>::operator=(
        TensorCopyFineToCoarseZerothOrder3D<T, nDim> const &rhs)
{
    dimDx = rhs.dimDx;
    dimDt = rhs.dimDt;
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T, int nDim>
void TensorCopyFineToCoarseZerothOrder3D<T, nDim>::process(
    Box3D fineDomain, TensorField3D<T, nDim> &fineField, TensorField3D<T, nDim> &coarseField)
{
    Dot3D posFine = fineField.getLocation();
    Dot3D posCoarse = coarseField.getLocation();

    plint stretch = util::twoToThePower(dxScale);
    Box3D coarseDomain(
        fineDomain.shift(posFine.x, posFine.y, posFine.z)
            .  // Convert to absolute fine coordinates.
        divideAndFitSmaller(stretch)
            .  // Rescale, but don't exceed original domain.
        shift(
            -posCoarse.x, -posCoarse.y, -posCoarse.z));  // Convert to relative coarse coordinates.
    PLB_ASSERT(contained(coarseDomain, coarseField.getBoundingBox()));

    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);

    plint fineX = (coarseDomain.x0 + posCoarse.x) * stretch - posFine.x;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX, fineX += stretch) {
        plint fineY = (coarseDomain.y0 + posCoarse.y) * stretch - posFine.y;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1;
             ++coarseY, fineY += stretch) {
            plint fineZ = (coarseDomain.z0 + posCoarse.z) * stretch - posFine.z;
            for (plint coarseZ = coarseDomain.z0; coarseZ <= coarseDomain.z1;
                 ++coarseZ, fineZ += stretch) {
                coarseField.get(coarseX, coarseY, coarseZ) =
                    fineField.get(fineX, fineY, fineZ) * scaleFactor;
            }
        }
    }
}

template <typename T, int nDim>
TensorCopyFineToCoarseZerothOrder3D<T, nDim> *TensorCopyFineToCoarseZerothOrder3D<T, nDim>::clone()
    const
{
    return new TensorCopyFineToCoarseZerothOrder3D<T, nDim>(*this);
}

/* ******************* LatticeCopyFineToCoarseZerothOrder3D ******************* */

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor>::LatticeCopyFineToCoarseZerothOrder3D(
    plint dxScale_, plint dtScale_) :
    dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be positive because the conversion goes from fine to coarse.
    PLB_ASSERT(dxScale >= 1);
}

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor>::LatticeCopyFineToCoarseZerothOrder3D(
    LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor> const &rhs) :
    dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor>
    &LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor>::operator=(
        LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor> const &rhs)
{
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor>::process(
    Box3D fineDomain, BlockLattice3D<T, Descriptor> &fineLattice,
    BlockLattice3D<T, Descriptor> &coarseLattice)
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

    Dot3D posFine = fineLattice.getLocation();
    Dot3D posCoarse = coarseLattice.getLocation();

    plint stretch = util::intTwoToThePower(dxScale);
    Box3D coarseDomain(
        fineDomain.shift(posFine.x, posFine.y, posFine.z)
            .  // Convert to absolute fine coordinates.
        divideAndFitSmaller(stretch)
            .  // Rescale, but don't exceed original domain.
        shift(
            -posCoarse.x, -posCoarse.y, -posCoarse.z));  // Convert to relative coarse coordinates.
    PLB_ASSERT(contained(coarseDomain, coarseLattice.getBoundingBox()));

    plint fineX = (coarseDomain.x0 + posCoarse.x) * stretch - posFine.x;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX, fineX += stretch) {
        plint fineY = (coarseDomain.y0 + posCoarse.y) * stretch - posFine.y;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1;
             ++coarseY, fineY += stretch) {
            plint fineZ = (coarseDomain.z0 + posCoarse.z) * stretch - posFine.z;
            for (plint coarseZ = coarseDomain.z0; coarseZ <= coarseDomain.z1;
                 ++coarseZ, fineZ += stretch) {
                // Fine --> Coarse
                Cell<T, Descriptor> const &fineCell = fineLattice.get(fineX, fineY, fineZ);
                Cell<T, Descriptor> &coarseCell = coarseLattice.get(coarseX, coarseY, coarseZ);
                T fineOmega = fineCell.getDynamics().getOmega();
                T fine_nu_cs2 = 1. / fineOmega - 0.5;
                T coarse_nu_cs2 = viscScale * fine_nu_cs2;
                T coarseOmega = 1. / (coarse_nu_cs2 + 0.5);
                coarseCell.getDynamics().setOmega(coarseOmega);

                T localFneqScale = fNeqScale * fineOmega / coarseOmega;
                plint numComponents = 1 + Descriptor<T>::d + Descriptor<T>::q
                                      + Descriptor<T>::ExternalField::numScalars;
                // TODO: For now external scalars are set to zero during conversion. They should be
                //   rescaled properly.
                std::vector<T> components(numComponents);
                std::fill(components.begin(), components.end(), 0.);
                plint order = 0;
                fineCell.getDynamics().decompose(fineCell, components, order);
                T rhoFine = Descriptor<T>::fullRho(components[0]);
                if (!fineCell.getDynamics().velIsJ()) {
                    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
                        components[1 + iDim] /= rhoFine;
                    }
                }
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
                // coarseCell.getDynamics().recompose(coarseCell, components, order);
                //  For recomposition, get the dynamics from the fine cell, because the fine/coarse
                //  dynamics may not have the same number of variables in the decomposed vector.
                fineCell.getDynamics().recompose(coarseCell, components, order);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor>
    *LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor>::clone() const
{
    return new LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor>(*this);
}

/* ******************* LatticeCopyFineToCoarseNoScale3D ******************* */

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseNoScale3D<T, Descriptor>::LatticeCopyFineToCoarseNoScale3D(plint dxScale_) :
    dxScale(dxScale_)
{
    // X scale must be positive because the conversion goes from fine to coarse.
    PLB_ASSERT(dxScale >= 1);
}

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseNoScale3D<T, Descriptor>::LatticeCopyFineToCoarseNoScale3D(
    LatticeCopyFineToCoarseNoScale3D<T, Descriptor> const &rhs) :
    dxScale(rhs.dxScale)
{ }

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseNoScale3D<T, Descriptor>
    &LatticeCopyFineToCoarseNoScale3D<T, Descriptor>::operator=(
        LatticeCopyFineToCoarseNoScale3D<T, Descriptor> const &rhs)
{
    dxScale = rhs.dxScale;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void LatticeCopyFineToCoarseNoScale3D<T, Descriptor>::process(
    Box3D fineDomain, BlockLattice3D<T, Descriptor> &fineLattice,
    BlockLattice3D<T, Descriptor> &coarseLattice)
{
    Dot3D posFine = fineLattice.getLocation();
    Dot3D posCoarse = coarseLattice.getLocation();

    plint stretch = util::twoToThePower(dxScale);
    Box3D coarseDomain(
        fineDomain.shift(posFine.x, posFine.y, posFine.z)
            .  // Convert to absolute fine coordinates.
        divideAndFitSmaller(stretch)
            .  // Rescale, but don't exceed original domain.
        shift(
            -posCoarse.x, -posCoarse.y, -posCoarse.z));  // Convert to relative coarse coordinates.
    PLB_ASSERT(contained(coarseDomain, coarseLattice.getBoundingBox()));

    plint fineX = (coarseDomain.x0 + posCoarse.x) * stretch - posFine.x;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX, fineX += stretch) {
        plint fineY = (coarseDomain.y0 + posCoarse.y) * stretch - posFine.y;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1;
             ++coarseY, fineY += stretch) {
            plint fineZ = (coarseDomain.z0 + posCoarse.z) * stretch - posFine.z;
            for (plint coarseZ = coarseDomain.z0; coarseZ <= coarseDomain.z1;
                 ++coarseZ, fineZ += stretch) {
                // Fine --> Coarse
                Cell<T, Descriptor> const &fineCell = fineLattice.get(fineX, fineY, fineZ);
                Cell<T, Descriptor> &coarseCell = coarseLattice.get(coarseX, coarseY, coarseZ);
                coarseCell.getDynamics().setOmega(fineCell.getDynamics().getOmega());
                coarseCell.attributeF(fineCell);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeCopyFineToCoarseNoScale3D<T, Descriptor>
    *LatticeCopyFineToCoarseNoScale3D<T, Descriptor>::clone() const
{
    return new LatticeCopyFineToCoarseNoScale3D<T, Descriptor>(*this);
}

/* ******************* CopyRescaleParticles3D ******************* */

template <typename T, template <typename U> class Descriptor>
CopyRescaleParticles3D<T, Descriptor>::CopyRescaleParticles3D(plint dxScale_, plint dtScale_) :
    dxScale(dxScale_), dtScale(dtScale_)
{ }

template <typename T, template <typename U> class Descriptor>
void CopyRescaleParticles3D<T, Descriptor>::processGenericBlocks(
    Box3D domain1, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ParticleField3D<T, Descriptor> &field1 =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[0]);
    ParticleField3D<T, Descriptor> &field2 =
        *dynamic_cast<ParticleField3D<T, Descriptor> *>(blocks[1]);

    Box3D domain2(field2.getBoundingBox());
    std::vector<Particle3D<T, Descriptor> *> particles;
    field1.findParticles(domain1, particles);
    for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
        Particle3D<T, Descriptor> *newparticle = particles[iParticle]->clone();
        newparticle->rescale(dxScale, dtScale);
        field2.addParticle(domain2, newparticle);
    }
}

template <typename T, template <typename U> class Descriptor>
CopyRescaleParticles3D<T, Descriptor> *CopyRescaleParticles3D<T, Descriptor>::clone() const
{
    return new CopyRescaleParticles3D<T, Descriptor>(*this);
}

/* ******************* ScalarInterpolateCoarseToFine3D ******************* */

template <typename T>
ScalarInterpolateCoarseToFine3D<T>::ScalarInterpolateCoarseToFine3D(
    plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_) :
    dimDx(dimDx_), dimDt(dimDt_), dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be negative because the conversion goes from coarse to fine.
    PLB_ASSERT(dxScale <= -1);
}

template <typename T>
ScalarInterpolateCoarseToFine3D<T>::ScalarInterpolateCoarseToFine3D(
    ScalarInterpolateCoarseToFine3D<T> const &rhs) :
    dimDx(rhs.dimDx), dimDt(rhs.dimDt), dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T>
ScalarInterpolateCoarseToFine3D<T> &ScalarInterpolateCoarseToFine3D<T>::operator=(
    ScalarInterpolateCoarseToFine3D<T> const &rhs)
{
    dimDx = rhs.dimDx;
    dimDt = rhs.dimDt;
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T>
void ScalarInterpolateCoarseToFine3D<T>::process(
    Box3D coarseDomain, ScalarField3D<T> &coarseField, ScalarField3D<T> &fineField)
{
    Dot3D posFine = fineField.getLocation();
    Dot3D posCoarse = coarseField.getLocation();

    plint stretch = util::twoToThePower(-dxScale);
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);

    std::vector<Dot3D> cellPos;
    std::vector<T> weights;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX) {
        plint fineX = (coarseX + posCoarse.x) * stretch - posFine.x;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1; ++coarseY) {
            plint fineY = (coarseY + posCoarse.y) * stretch - posFine.y;
            for (plint coarseZ = coarseDomain.z0; coarseZ <= coarseDomain.z1; ++coarseZ) {
                plint fineZ = (coarseZ + posCoarse.z) * stretch - posFine.z;
                for (plint iX = fineX - stretch / 2 + 1; iX <= fineX + stretch / 2; ++iX) {
                    for (plint iY = fineY - stretch / 2 + 1; iY <= fineY + stretch / 2; ++iY) {
                        for (plint iZ = fineZ - stretch / 2 + 1; iZ <= fineZ + stretch / 2; ++iZ) {
                            if (contained(iX, iY, iZ, fineField.getBoundingBox())) {
                                T coarseIx = (T)(iX + posFine.x) / (T)stretch;
                                T coarseIy = (T)(iY + posFine.y) / (T)stretch;
                                T coarseIz = (T)(iZ + posFine.z) / (T)stretch;
                                linearInterpolationCoefficients<T>(
                                    coarseField, Array<T, 3>(coarseIx, coarseIy, coarseIz), cellPos,
                                    weights);
                                fineField.get(iX, iY, iZ) = 0.;
                                for (plint i = 0; i < 8; ++i) {
                                    fineField.get(iX, iY, iZ) +=
                                        weights[i]
                                        * coarseField.get(cellPos[i].x, cellPos[i].y, cellPos[i].z);
                                }
                                fineField.get(iX, iY, iZ) *= scaleFactor;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename T>
ScalarInterpolateCoarseToFine3D<T> *ScalarInterpolateCoarseToFine3D<T>::clone() const
{
    return new ScalarInterpolateCoarseToFine3D<T>(*this);
}

/* ******************* TensorInterpolateCoarseToFine3D ******************* */

template <typename T, int nDim>
TensorInterpolateCoarseToFine3D<T, nDim>::TensorInterpolateCoarseToFine3D(
    plint dimDx_, plint dimDt_, plint dxScale_, plint dtScale_) :
    dimDx(dimDx_), dimDt(dimDt_), dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be negative because the conversion goes from coarse to fine.
    PLB_ASSERT(dxScale <= -1);
}

template <typename T, int nDim>
TensorInterpolateCoarseToFine3D<T, nDim>::TensorInterpolateCoarseToFine3D(
    TensorInterpolateCoarseToFine3D<T, nDim> const &rhs) :
    dimDx(rhs.dimDx), dimDt(rhs.dimDt), dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T, int nDim>
TensorInterpolateCoarseToFine3D<T, nDim> &TensorInterpolateCoarseToFine3D<T, nDim>::operator=(
    TensorInterpolateCoarseToFine3D<T, nDim> const &rhs)
{
    dimDx = rhs.dimDx;
    dimDt = rhs.dimDt;
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T, int nDim>
void TensorInterpolateCoarseToFine3D<T, nDim>::process(
    Box3D coarseDomain, TensorField3D<T, nDim> &coarseField, TensorField3D<T, nDim> &fineField)
{
    Dot3D posFine = fineField.getLocation();
    Dot3D posCoarse = coarseField.getLocation();

    plint stretch = util::twoToThePower(-dxScale);
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);

    std::vector<Dot3D> cellPos;
    std::vector<T> weights;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX) {
        plint fineX = (coarseX + posCoarse.x) * stretch - posFine.x;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1; ++coarseY) {
            plint fineY = (coarseY + posCoarse.y) * stretch - posFine.y;
            for (plint coarseZ = coarseDomain.z0; coarseZ <= coarseDomain.z1; ++coarseZ) {
                plint fineZ = (coarseZ + posCoarse.z) * stretch - posFine.z;
                for (plint iX = fineX - stretch / 2 + 1; iX <= fineX + stretch / 2; ++iX) {
                    for (plint iY = fineY - stretch / 2 + 1; iY <= fineY + stretch / 2; ++iY) {
                        for (plint iZ = fineZ - stretch / 2 + 1; iZ <= fineZ + stretch / 2; ++iZ) {
                            if (contained(iX, iY, iZ, fineField.getBoundingBox())) {
                                T coarseIx = (T)(iX + posFine.x) / (T)stretch;
                                T coarseIy = (T)(iY + posFine.y) / (T)stretch;
                                T coarseIz = (T)(iZ + posFine.z) / (T)stretch;
                                linearInterpolationCoefficients(
                                    coarseField, Array<T, 3>(coarseIx, coarseIy, coarseIz), cellPos,
                                    weights);
                                fineField.get(iX, iY, iZ).resetToZero();
                                for (plint i = 0; i < 8; ++i) {
                                    fineField.get(iX, iY, iZ) +=
                                        weights[i]
                                        * coarseField.get(cellPos[i].x, cellPos[i].y, cellPos[i].z);
                                }
                                fineField.get(iX, iY, iZ) *= scaleFactor;
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename T, int nDim>
TensorInterpolateCoarseToFine3D<T, nDim> *TensorInterpolateCoarseToFine3D<T, nDim>::clone() const
{
    return new TensorInterpolateCoarseToFine3D<T, nDim>(*this);
}

/* ******************* LatticeInterpolateCoarseToFine3D ******************* */

template <typename T, template <typename U> class Descriptor>
LatticeInterpolateCoarseToFine3D<T, Descriptor>::LatticeInterpolateCoarseToFine3D(
    plint dxScale_, plint dtScale_) :
    dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be negative because the conversion goes from coarse to fine.
    PLB_ASSERT(dxScale <= -1);
}

template <typename T, template <typename U> class Descriptor>
LatticeInterpolateCoarseToFine3D<T, Descriptor>::LatticeInterpolateCoarseToFine3D(
    LatticeInterpolateCoarseToFine3D<T, Descriptor> const &rhs) :
    dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T, template <typename U> class Descriptor>
LatticeInterpolateCoarseToFine3D<T, Descriptor>
    &LatticeInterpolateCoarseToFine3D<T, Descriptor>::operator=(
        LatticeInterpolateCoarseToFine3D<T, Descriptor> const &rhs)
{
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void LatticeInterpolateCoarseToFine3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
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

    Dot3D posFine = fineLattice.getLocation();
    Dot3D posCoarse = coarseLattice.getLocation();

    plint stretch = util::intTwoToThePower(-dxScale);

    std::vector<Dot3D> cellPos;
    std::vector<T> weights;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX) {
        plint fineX = (coarseX + posCoarse.x) * stretch - posFine.x;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1; ++coarseY) {
            plint fineY = (coarseY + posCoarse.y) * stretch - posFine.y;
            for (plint coarseZ = coarseDomain.z0; coarseZ <= coarseDomain.z1; ++coarseZ) {
                plint fineZ = (coarseZ + posCoarse.z) * stretch - posFine.z;
                for (plint iX = fineX - stretch / 2 + 1; iX <= fineX + stretch / 2; ++iX) {
                    for (plint iY = fineY - stretch / 2 + 1; iY <= fineY + stretch / 2; ++iY) {
                        for (plint iZ = fineZ - stretch / 2 + 1; iZ <= fineZ + stretch / 2; ++iZ) {
                            if (contained(iX, iY, iZ, fineLattice.getBoundingBox())) {
                                T coarseIx = (T)(iX + posFine.x) / (T)stretch;
                                T coarseIy = (T)(iY + posFine.y) / (T)stretch;
                                T coarseIz = (T)(iZ + posFine.z) / (T)stretch;
                                linearInterpolationCoefficients(
                                    coarseLattice, Array<T, 3>(coarseIx, coarseIy, coarseIz),
                                    cellPos, weights);
                                plint order = 0;
                                plint numComponents = 1 + Descriptor<T>::d + Descriptor<T>::q
                                                      + Descriptor<T>::ExternalField::numScalars;
                                // TODO: For now external scalars are set to zero during conversion.
                                // They should be
                                //   rescaled properly.
                                std::vector<T> components(numComponents),
                                    oldComponents(numComponents);
                                std::fill(components.begin(), components.end(), 0.);
                                // T coarse_nu_cs2 = 0.;
                                components[0] =
                                    1 - Descriptor<T>::SkordosFactor();  // Density offset.
                                for (plint i = 0; i < 8; ++i) {
                                    Cell<T, Descriptor> const &coarseCell =
                                        coarseLattice.get(cellPos[i].x, cellPos[i].y, cellPos[i].z);
                                    // T coarseOmega = coarseCell.getDynamics().getOmega();
                                    // T tmp_coarse_nu_cs2 = 1./coarseOmega - 0.5;
                                    // coarse_nu_cs2 += weights[i]*tmp_coarse_nu_cs2;
                                    coarseCell.getDynamics().decompose(
                                        coarseCell, oldComponents, order);
                                    // 1. rhoBar. Keep the constant term (1), and rescale the
                                    // remaining term like a pressure.
                                    T pres_cs2 =
                                        oldComponents[0] - (1 - Descriptor<T>::SkordosFactor());
                                    components[0] += weights[i] * pres_cs2 * presScale;
                                    // 2. j.
                                    T velMultFactor = weights[i] * velScale;
                                    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
                                        components[1 + iDim] +=
                                            velMultFactor * oldComponents[1 + iDim];
                                    }
                                    // 3. Off-equilibrium populations.
                                    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                                        components[1 + Descriptor<T>::d + iPop] +=
                                            weights[i] * fNeqScale
                                            * oldComponents[1 + Descriptor<T>::d + iPop];
                                    }
                                }
                                Cell<T, Descriptor> &fineCell = fineLattice.get(iX, iY, iZ);
                                // Rescaling of the relaxation time has been deactivated for now,
                                // because it produces numerical instability when a checkpoint file
                                // is read and scaled. Needs to be fixed at some point.
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
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeInterpolateCoarseToFine3D<T, Descriptor>
    *LatticeInterpolateCoarseToFine3D<T, Descriptor>::clone() const
{
    return new LatticeInterpolateCoarseToFine3D<T, Descriptor>(*this);
}

/* ******************* LatticeLinearInterpolateCoarseToFine3D ******************* */

template <typename T, template <typename U> class Descriptor>
LatticeLinearInterpolateCoarseToFine3D<T, Descriptor>::LatticeLinearInterpolateCoarseToFine3D(
    plint dxScale_, plint dtScale_) :
    dxScale(dxScale_), dtScale(dtScale_)
{
    // X scale must be negative because the conversion goes from coarse to fine.
    PLB_ASSERT(dxScale <= -1);
}

template <typename T, template <typename U> class Descriptor>
LatticeLinearInterpolateCoarseToFine3D<T, Descriptor>::LatticeLinearInterpolateCoarseToFine3D(
    LatticeLinearInterpolateCoarseToFine3D<T, Descriptor> const &rhs) :
    dxScale(rhs.dxScale), dtScale(rhs.dtScale)
{ }

template <typename T, template <typename U> class Descriptor>
LatticeLinearInterpolateCoarseToFine3D<T, Descriptor>
    &LatticeLinearInterpolateCoarseToFine3D<T, Descriptor>::operator=(
        LatticeLinearInterpolateCoarseToFine3D<T, Descriptor> const &rhs)
{
    dxScale = rhs.dxScale;
    dtScale = rhs.dtScale;
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void LatticeLinearInterpolateCoarseToFine3D<T, Descriptor>::process(
    Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
    BlockLattice3D<T, Descriptor> &fineLattice)
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

    Dot3D posFine = fineLattice.getLocation();
    Dot3D posCoarse = coarseLattice.getLocation();

    plint stretch = util::intTwoToThePower(-dxScale);

    std::vector<Dot3D> cellPos;
    std::vector<T> weights;
    for (plint coarseX = coarseDomain.x0; coarseX <= coarseDomain.x1; ++coarseX) {
        plint fineX = (coarseX + posCoarse.x) * stretch - posFine.x;
        for (plint coarseY = coarseDomain.y0; coarseY <= coarseDomain.y1; ++coarseY) {
            plint fineY = (coarseY + posCoarse.y) * stretch - posFine.y;
            for (plint coarseZ = coarseDomain.z0; coarseZ <= coarseDomain.z1; ++coarseZ) {
                plint fineZ = (coarseZ + posCoarse.z) * stretch - posFine.z;
                for (plint iX = fineX - stretch / 2 + 1; iX <= fineX + stretch / 2; ++iX) {
                    for (plint iY = fineY - stretch / 2 + 1; iY <= fineY + stretch / 2; ++iY) {
                        for (plint iZ = fineZ - stretch / 2 + 1; iZ <= fineZ + stretch / 2; ++iZ) {
                            if (contained(iX, iY, iZ, fineLattice.getBoundingBox())) {
                                T coarseIx = (T)(iX + posFine.x) / (T)stretch;
                                T coarseIy = (T)(iY + posFine.y) / (T)stretch;
                                T coarseIz = (T)(iZ + posFine.z) / (T)stretch;
                                linearInterpolationCoefficients(
                                    coarseLattice, Array<T, 3>(coarseIx, coarseIy, coarseIz),
                                    cellPos, weights);
                                plint order = 0;
                                plint numComponents = 1 + Descriptor<T>::d + Descriptor<T>::q
                                                      + Descriptor<T>::ExternalField::numScalars;
                                // TODO: For now external scalars are set to zero during conversion.
                                // They should be
                                //   rescaled properly.
                                std::vector<T> components(numComponents),
                                    oldComponents(numComponents);
                                std::fill(components.begin(), components.end(), 0.);
                                T coarse_nu_cs2 = 0.;
                                components[0] =
                                    1 - Descriptor<T>::SkordosFactor();  // Density offset.
                                for (plint i = 0; i < 8; ++i) {
                                    Cell<T, Descriptor> const &coarseCell =
                                        coarseLattice.get(cellPos[i].x, cellPos[i].y, cellPos[i].z);
                                    T coarseOmega = coarseCell.getDynamics().getOmega();
                                    T tmp_coarse_nu_cs2 = 1. / coarseOmega - 0.5;
                                    coarse_nu_cs2 += weights[i] * tmp_coarse_nu_cs2;
                                    coarseCell.getDynamics().decompose(
                                        coarseCell, oldComponents, order);
                                    // 1. rhoBar. Keep the constant term (1), and rescale the
                                    // remaining term like a pressure.
                                    T pres_cs2 =
                                        oldComponents[0] - (1 - Descriptor<T>::SkordosFactor());
                                    components[0] += weights[i] * pres_cs2 * presScale;
                                    // 2. j.
                                    T velMultFactor = weights[i] * velScale;
                                    for (plint iDim = 0; iDim < Descriptor<T>::d; ++iDim) {
                                        components[1 + iDim] +=
                                            velMultFactor * oldComponents[1 + iDim];
                                    }
                                    // 3. Off-equilibrium populations.
                                    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                                        components[1 + Descriptor<T>::d + iPop] +=
                                            weights[i] * fNeqScale
                                            * oldComponents[1 + Descriptor<T>::d + iPop];
                                    }
                                }
                                Cell<T, Descriptor> &fineCell = fineLattice.get(iX, iY, iZ);
                                // Rescaling of the relaxation time has been deactivated for now,
                                // because it produces numerical instability when a checkpoint file
                                // is read and scaled. Needs to be fixed at some point.

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
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeLinearInterpolateCoarseToFine3D<T, Descriptor>
    *LatticeLinearInterpolateCoarseToFine3D<T, Descriptor>::clone() const
{
    return new LatticeLinearInterpolateCoarseToFine3D<T, Descriptor>(*this);
}

/* ******************* Free functions ******************* */

inline Box3D scaleBox(Box3D box, plint nLevel)
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

inline MultiBlockManagement3D scaleMultiBlockManagement(
    MultiBlockManagement3D const &multiBlockManagement, plint nLevel)
{
    SparseBlockStructure3D const &sparseBlock = multiBlockManagement.getSparseBlockStructure();

    Box3D rescaledBoundingBox = scaleBox(sparseBlock.getBoundingBox(), nLevel);
    SparseBlockStructure3D scaledSparseBlock(rescaledBoundingBox);
    plint stretch = 1;
    if (nLevel > 0) {
        stretch = util::intTwoToThePower(nLevel);
    }

    std::map<plint, Box3D>::const_iterator it = sparseBlock.getBulks().begin();
    for (; it != sparseBlock.getBulks().end(); ++it) {
        Box3D scaledBulk = scaleBox(it->second, nLevel);
        if (stretch > 1) {
            scaledBulk.x0 -= stretch / 2 - 1;
            scaledBulk.y0 -= stretch / 2 - 1;
            scaledBulk.z0 -= stretch / 2 - 1;
            scaledBulk.x1 += stretch / 2;
            scaledBulk.y1 += stretch / 2;
            scaledBulk.z1 += stretch / 2;
        }
        Box3D scaledUniqueBulk = scaledBulk;  // TODO: compute unique bulk properly.
        plint blockId = it->first;
        scaledSparseBlock.addBlock(scaledBulk, scaledUniqueBulk, blockId);
    }
    return MultiBlockManagement3D(
        scaledSparseBlock, multiBlockManagement.getThreadAttribution().clone(),
        multiBlockManagement.getEnvelopeWidth(),
        multiBlockManagement.getRefinementLevel() + nLevel);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > coarsen(
    MultiScalarField3D<T> &fineField, plint dimDx, plint dimDt, plint dxScale, plint dtScale)
{
    PLB_PRECONDITION(dxScale >= 1);
    // Relative level is negative when going from fine to coarse.
    plint relativeLevel = -dxScale;
    MultiBlockManagement3D management =
        scaleMultiBlockManagement(fineField.getMultiBlockManagement(), relativeLevel);
    MultiScalarField3D<T> *result = new MultiScalarField3D<T>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>());
    applyProcessingFunctional(
        new ScalarCopyFineToCoarseZerothOrder3D<T>(dimDx, dimDt, dxScale, dtScale),
        fineField.getBoundingBox(), fineField, *result);

    return std::unique_ptr<MultiScalarField3D<T> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > coarsen(
    MultiTensorField3D<T, nDim> &fineField, plint dimDx, plint dimDt, plint dxScale, plint dtScale)
{
    PLB_PRECONDITION(dxScale >= 1);
    // Relative level is negative when going from fine to coarse.
    plint relativeLevel = -dxScale;
    MultiBlockManagement3D management =
        scaleMultiBlockManagement(fineField.getMultiBlockManagement(), relativeLevel);
    MultiTensorField3D<T, nDim> *result = new MultiTensorField3D<T, nDim>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>());
    applyProcessingFunctional(
        new TensorCopyFineToCoarseZerothOrder3D<T, nDim>(dimDx, dimDt, dxScale, dtScale),
        fineField.getBoundingBox(), fineField, *result);

    return std::unique_ptr<MultiTensorField3D<T, nDim> >(result);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > refine(
    MultiScalarField3D<T> &coarseField, plint dimDx, plint dimDt, plint dxScale, plint dtScale)
{
    PLB_PRECONDITION(dxScale <= -1);
    // Relative level is positive when going from coarse to fine.
    plint relativeLevel = -dxScale;
    MultiBlockManagement3D management =
        scaleMultiBlockManagement(coarseField.getMultiBlockManagement(), relativeLevel);
    MultiScalarField3D<T> *result = new MultiScalarField3D<T>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiScalarAccess<T>());
    applyProcessingFunctional(
        new ScalarInterpolateCoarseToFine3D<T>(dimDx, dimDt, dxScale, dtScale),
        coarseField.getBoundingBox(), coarseField, *result);

    return std::unique_ptr<MultiScalarField3D<T> >(result);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > refine(
    MultiTensorField3D<T, nDim> &coarseField, plint dimDx, plint dimDt, plint dxScale,
    plint dtScale)
{
    PLB_PRECONDITION(dxScale <= -1);
    // Relative level is positive when going from coarse to fine.
    plint relativeLevel = -dxScale;
    MultiBlockManagement3D management =
        scaleMultiBlockManagement(coarseField.getMultiBlockManagement(), relativeLevel);
    MultiTensorField3D<T, nDim> *result = new MultiTensorField3D<T, nDim>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>());
    applyProcessingFunctional(
        new TensorInterpolateCoarseToFine3D<T, nDim>(dimDx, dimDt, dxScale, dtScale),
        coarseField.getBoundingBox(), coarseField, *result);

    return std::unique_ptr<MultiTensorField3D<T, nDim> >(result);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > coarsen(
    MultiBlockLattice3D<T, Descriptor> &fineLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics)
{
    PLB_PRECONDITION(dxScale >= 1);
    // Relative level is negative when going from fine to coarse.
    plint relativeLevel = -dxScale;
    MultiBlockManagement3D management =
        scaleMultiBlockManagement(fineLattice.getMultiBlockManagement(), relativeLevel);
    MultiBlockLattice3D<T, Descriptor> *result = new MultiBlockLattice3D<T, Descriptor>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
        backgroundDynamics->clone());
    defineDynamics(*result, result->getBoundingBox(), backgroundDynamics);
    applyProcessingFunctional(
        new LatticeCopyFineToCoarseZerothOrder3D<T, Descriptor>(dxScale, dtScale),
        fineLattice.getBoundingBox(), fineLattice, *result);

    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(result);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > coarsenNoScale(
    MultiBlockLattice3D<T, Descriptor> &fineLattice, plint dxScale,
    Dynamics<T, Descriptor> *backgroundDynamics)
{
    PLB_PRECONDITION(dxScale >= 1);
    // Relative level is negative when going from fine to coarse.
    plint relativeLevel = -dxScale;
    MultiBlockManagement3D management =
        scaleMultiBlockManagement(fineLattice.getMultiBlockManagement(), relativeLevel);
    MultiBlockLattice3D<T, Descriptor> *result = new MultiBlockLattice3D<T, Descriptor>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
        backgroundDynamics->clone());
    defineDynamics(*result, result->getBoundingBox(), backgroundDynamics);
    applyProcessingFunctional(
        new LatticeCopyFineToCoarseNoScale3D<T, Descriptor>(dxScale), fineLattice.getBoundingBox(),
        fineLattice, *result);

    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(result);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > refine(
    MultiBlockLattice3D<T, Descriptor> &coarseLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics)
{
    PLB_PRECONDITION(dxScale <= -1);
    // Relative level is positive when going from coarse to fine.
    plint relativeLevel = -dxScale;
    MultiBlockManagement3D management =
        scaleMultiBlockManagement(coarseLattice.getMultiBlockManagement(), relativeLevel);
    MultiBlockLattice3D<T, Descriptor> *result = new MultiBlockLattice3D<T, Descriptor>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
        backgroundDynamics->clone());
    defineDynamics(*result, result->getBoundingBox(), backgroundDynamics);
    applyProcessingFunctional(
        new LatticeInterpolateCoarseToFine3D<T, Descriptor>(dxScale, dtScale),
        coarseLattice.getBoundingBox(), coarseLattice, *result);

    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(result);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > linearRefine(
    MultiBlockLattice3D<T, Descriptor> &coarseLattice, plint dxScale, plint dtScale,
    Dynamics<T, Descriptor> *backgroundDynamics)
{
    PLB_PRECONDITION(dxScale <= -1);
    // Relative level is positive when going from coarse to fine.
    plint relativeLevel = -dxScale;
    MultiBlockManagement3D management =
        scaleMultiBlockManagement(coarseLattice.getMultiBlockManagement(), relativeLevel);
    MultiBlockLattice3D<T, Descriptor> *result = new MultiBlockLattice3D<T, Descriptor>(
        management, defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
        backgroundDynamics->clone());
    defineDynamics(*result, result->getBoundingBox(), backgroundDynamics);
    applyProcessingFunctional(
        new LatticeLinearInterpolateCoarseToFine3D<T, Descriptor>(dxScale, dtScale),
        coarseLattice.getBoundingBox(), coarseLattice, *result);

    return std::unique_ptr<MultiBlockLattice3D<T, Descriptor> >(result);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > > refine(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &coarseField, plint dxScale,
    plint dtScale)
{
    PLB_PRECONDITION(dxScale <= -1);
    // Relative level is positive when going from coarse to fine.
    plint relativeLevel = -dxScale;
    MultiBlockManagement3D management =
        scaleMultiBlockManagement(coarseField.getMultiBlockManagement(), relativeLevel);
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > *result =
        new MultiParticleField3D<DenseParticleField3D<T, Descriptor> >(
            management, defaultMultiBlockPolicy3D().getCombinedStatistics());
    std::vector<MultiBlock3D *> args;
    args.push_back(&coarseField);
    args.push_back(result);
    applyProcessingFunctional(
        new CopyRescaleParticles3D<T, Descriptor>(dxScale, dtScale), coarseField.getBoundingBox(),
        args);

    return std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > >(result);
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > > coarsen(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &fineField, plint dxScale,
    plint dtScale)
{
    PLB_PRECONDITION(dxScale >= 1);
    // Relative level is negative when going from fine to coarse.
    plint relativeLevel = -dxScale;
    MultiBlockManagement3D management =
        scaleMultiBlockManagement(fineField.getMultiBlockManagement(), relativeLevel);
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > *result =
        new MultiParticleField3D<DenseParticleField3D<T, Descriptor> >(
            management, defaultMultiBlockPolicy3D().getCombinedStatistics());

    std::vector<MultiBlock3D *> args;
    args.push_back(&fineField);
    args.push_back(result);
    applyProcessingFunctional(
        new CopyRescaleParticles3D<T, Descriptor>(dxScale, dtScale), fineField.getBoundingBox(),
        args);

    return std::unique_ptr<MultiParticleField3D<DenseParticleField3D<T, Descriptor> > >(result);
}

}  // namespace plb

#endif  // GRID_CONVERSION_3D_HH
