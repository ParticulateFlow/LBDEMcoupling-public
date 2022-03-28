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
 * Multi scale manager -- generic implementation.
 */

#include "multiGrid/multiScale.h"

#include "multiBlock/multiBlockManagement2D.h"
#include "multiBlock/multiBlockManagement3D.h"

namespace plb {

///////////////// Class PowerTwoMultiScaleManager //////////////////////////

plint PowerTwoMultiScaleManager::twoToTheLevel(plint nLevel)
{
    plint powerOfTwo = 1;
    for (plint iL = 0; iL < nLevel; ++iL) {
        powerOfTwo *= 2;
    }
    return powerOfTwo;
}

Box2D PowerTwoMultiScaleManager::scaleBox(Box2D box, plint nLevel) const
{
    if (nLevel > 0) {
        return box.multiply(twoToTheLevel(nLevel));
    } else if (nLevel < 0) {
        // If the fine-grid box does not fit with the coarse-grid box,
        //   let's shrink the coarse-grid box by one fine cell. This
        //   makes sure we'll never exceed the fine-grid box
        return box.divideAndFitSmaller(twoToTheLevel(-nLevel));
    } else {
        return box;
    }
}

Box3D PowerTwoMultiScaleManager::scaleBox(Box3D box, plint nLevel) const
{
    if (nLevel > 0) {
        return box.multiply(twoToTheLevel(nLevel));
    } else if (nLevel < 0) {
        // If the fine-grid box does not fit with the coarse-grid box,
        //   let's shrink the coarse-grid box by one fine cell. This
        //   makes sure we'll never exceed the fine-grid box
        return box.divideAndFitSmaller(twoToTheLevel(-nLevel));
    } else {
        return box;
    }
}

MultiBlockManagement2D PowerTwoMultiScaleManager::scaleMultiBlockManagement(
    MultiBlockManagement2D const &multiBlockManagement, plint nLevel) const
{
    SparseBlockStructure2D const &sparseBlock = multiBlockManagement.getSparseBlockStructure();

    Box2D rescaledBoundingBox = scaleBox(sparseBlock.getBoundingBox(), nLevel);
    SparseBlockStructure2D scaledSparseBlock(rescaledBoundingBox);

    std::map<plint, Box2D>::const_iterator it = sparseBlock.getBulks().begin();
    for (; it != sparseBlock.getBulks().end(); ++it) {
        Box2D scaledBulk = scaleBox(it->second, nLevel);
        Box2D scaledUniqueBulk = scaledBulk;  // TODO: compute unique bulk properly.
        plint blockId = it->first;
        scaledSparseBlock.addBlock(scaledBulk, scaledUniqueBulk, blockId);
    }
    return MultiBlockManagement2D(
        scaledSparseBlock, multiBlockManagement.getThreadAttribution().clone(),
        multiBlockManagement.getEnvelopeWidth(),
        multiBlockManagement.getRefinementLevel() + nLevel);
}

MultiBlockManagement3D PowerTwoMultiScaleManager::scaleMultiBlockManagement(
    MultiBlockManagement3D const &multiBlockManagement, plint nLevel) const
{
    SparseBlockStructure3D const &sparseBlock = multiBlockManagement.getSparseBlockStructure();

    Box3D rescaledBoundingBox = scaleBox(sparseBlock.getBoundingBox(), nLevel);
    SparseBlockStructure3D scaledSparseBlock(rescaledBoundingBox);

    std::map<plint, Box3D>::const_iterator it = sparseBlock.getBulks().begin();
    for (; it != sparseBlock.getBulks().end(); ++it) {
        Box3D scaledBulk = scaleBox(it->second, nLevel);
        Box3D scaledUniqueBulk = scaledBulk;  // TODO: compute unique bulk properly.
        plint blockId = it->first;
        scaledSparseBlock.addBlock(scaledBulk, scaledUniqueBulk, blockId);
    }

    return MultiBlockManagement3D(
        scaledSparseBlock, multiBlockManagement.getThreadAttribution().clone(),
        multiBlockManagement.getEnvelopeWidth(),
        multiBlockManagement.getRefinementLevel() + nLevel);
}

///////////////// Class ConvectiveMultiScaleManager //////////////////////////

void ConvectiveMultiScaleManager::scaleVelocity(Array<double, 2> &u, plint nLevel) const
{
    // Velocity is scale-invariant in convective scaling.
}

void ConvectiveMultiScaleManager::scaleVelocity(Array<double, 3> &u, plint nLevel) const
{
    // Velocity is scale-invariant in convective scaling.
}

double ConvectiveMultiScaleManager::scaleDeltaX(plint nLevel) const
{
    if (nLevel > 0) {
        return 1. / (double)this->twoToTheLevel(nLevel);
    } else if (nLevel < 0) {
        return (double)this->twoToTheLevel(-nLevel);
    } else {
        return 1.;
    }
}

double ConvectiveMultiScaleManager::scaleDeltaT(plint nLevel) const
{
    // dt scale is equal to dx scale by definition in convective scaling.
    return scaleDeltaX(nLevel);
}

ConvectiveMultiScaleManager *ConvectiveMultiScaleManager::clone() const
{
    return new ConvectiveMultiScaleManager(*this);
}

namespace global {

DefaultMultiScaleManager::DefaultMultiScaleManager()
{
    defaultManager = new ConvectiveMultiScaleManager();
}

DefaultMultiScaleManager::~DefaultMultiScaleManager()
{
    delete defaultManager;
}

void DefaultMultiScaleManager::set(MultiScaleManager *newManager)
{
    delete defaultManager;
    defaultManager = newManager;
}

MultiScaleManager const &DefaultMultiScaleManager::get() const
{
    return *defaultManager;
}

DefaultMultiScaleManager &accessDefaultMultiScaleManager()
{
    static DefaultMultiScaleManager defaultMultiScaleManager;
    return defaultMultiScaleManager;
}

MultiScaleManager const &getDefaultMultiScaleManager()
{
    return accessDefaultMultiScaleManager().get();
}

void setDefaultMultiScaleManager(MultiScaleManager *newMultiScaleManager)
{
    accessDefaultMultiScaleManager().set(newMultiScaleManager);
}

}  // namespace global

}  // namespace plb
