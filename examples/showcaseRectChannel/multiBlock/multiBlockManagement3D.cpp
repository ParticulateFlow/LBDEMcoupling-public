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
 * Geometry specifications for 3D multiblocks -- implementation.
 */

#include "multiBlock/multiBlockManagement3D.h"

#include <algorithm>

#include "core/plbDebug.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/staticRepartitions3D.h"

namespace plb {

MultiBlockManagement3D::MultiBlockManagement3D(
    SparseBlockStructure3D const &sparseBlock_, ThreadAttribution *threadAttribution_,
    plint envelopeWidth_, plint refinementLevel_) :
    envelopeWidth(envelopeWidth_),
    sparseBlock(sparseBlock_),
    threadAttribution(threadAttribution_),
    localInfo(sparseBlock, getThreadAttribution(), envelopeWidth),
    refinementLevel(refinementLevel_)
{ }

MultiBlockManagement3D::MultiBlockManagement3D(MultiBlockManagement3D const &rhs) :
    envelopeWidth(rhs.envelopeWidth),
    sparseBlock(rhs.sparseBlock),
    threadAttribution(rhs.threadAttribution->clone()),
    localInfo(rhs.localInfo),
    refinementLevel(rhs.refinementLevel)
{ }

MultiBlockManagement3D &MultiBlockManagement3D::operator=(MultiBlockManagement3D const &rhs)
{
    MultiBlockManagement3D newBlockManagement(rhs);
    newBlockManagement.swap(*this);
    return *this;
}

void MultiBlockManagement3D::swap(MultiBlockManagement3D &rhs)
{
    std::swap(envelopeWidth, rhs.envelopeWidth);
    sparseBlock.swap(rhs.sparseBlock);
    std::swap(threadAttribution, rhs.threadAttribution);
    localInfo.swap(rhs.localInfo);
    std::swap(refinementLevel, rhs.refinementLevel);
}

MultiBlockManagement3D::~MultiBlockManagement3D()
{
    delete threadAttribution;
}

plint MultiBlockManagement3D::getEnvelopeWidth() const
{
    return envelopeWidth;
}

Box3D MultiBlockManagement3D::getBoundingBox() const
{
    return sparseBlock.getBoundingBox();
}

Box3D MultiBlockManagement3D::getBulk(plint blockId) const
{
    Box3D bulk;
#ifdef PLB_DEBUG
    bool ok =
#endif
        sparseBlock.getBulk(blockId, bulk);
    PLB_ASSERT(ok);
    return bulk;
}

Box3D MultiBlockManagement3D::getUniqueBulk(plint blockId) const
{
    Box3D uniqueBulk;
#ifdef PLB_DEBUG
    bool ok =
#endif
        sparseBlock.getUniqueBulk(blockId, uniqueBulk);
    PLB_ASSERT(ok);
    return uniqueBulk;
}

Box3D MultiBlockManagement3D::getEnvelope(plint blockId) const
{
    Box3D bulk;
#ifdef PLB_DEBUG
    bool ok =
#endif
        sparseBlock.getBulk(blockId, bulk);
    PLB_ASSERT(ok);
    return SmartBulk3D(sparseBlock, envelopeWidth, bulk).computeEnvelope();
}

SparseBlockStructure3D const &MultiBlockManagement3D::getSparseBlockStructure() const
{
    return sparseBlock;
}

LocalMultiBlockInfo3D const &MultiBlockManagement3D::getLocalInfo() const
{
    return localInfo;
}

ThreadAttribution const &MultiBlockManagement3D::getThreadAttribution() const
{
    return *threadAttribution;
}

void MultiBlockManagement3D::setCoProcessors(std::map<plint, int> const &coProcessors)
{
    threadAttribution->setCoProcessors(coProcessors);
}

bool MultiBlockManagement3D::findInLocalBulk(
    plint iX, plint iY, plint iZ, plint &foundId, plint &localX, plint &localY, plint &localZ) const
{
    foundId = sparseBlock.locate(iX, iY, iZ);
    SmartBulk3D bulk(sparseBlock, envelopeWidth, foundId);
    localX = bulk.toLocalX(iX);
    localY = bulk.toLocalY(iY);
    localZ = bulk.toLocalZ(iZ);
    return foundId >= 0;
}

bool MultiBlockManagement3D::findAllLocalRepresentations(
    plint iX, plint iY, plint iZ, std::vector<plint> &foundId, std::vector<plint> &foundX,
    std::vector<plint> &foundY, std::vector<plint> &foundZ) const
{
    bool hasBulkCell = false;
    // First, search in all blocks which are local to the current processor,
    // including in the envelopes.
    for (pluint iBlock = 0; iBlock < localInfo.getBlocks().size(); ++iBlock) {
        plint blockId = localInfo.getBlocks()[iBlock];
        SmartBulk3D bulk(sparseBlock, envelopeWidth, blockId);
        if (contained(iX, iY, iZ, bulk.computeEnvelope())) {
            if (contained(iX, iY, iZ, bulk.getBulk())) {
                hasBulkCell = true;
                foundId.insert(foundId.begin(), blockId);
                foundX.insert(foundX.begin(), bulk.toLocalX(iX));
                foundY.insert(foundY.begin(), bulk.toLocalY(iY));
                foundZ.insert(foundZ.begin(), bulk.toLocalZ(iZ));
            } else {
                foundId.push_back(blockId);
                foundX.push_back(bulk.toLocalX(iX));
                foundY.push_back(bulk.toLocalY(iY));
                foundZ.push_back(bulk.toLocalZ(iZ));
            }
        }
    }
    // Here's a subtlety: with periodic boundary conditions, one may need to take into
    //   account a cell which is not inside the boundingBox, because it's at the opposite
    //   boundary. Therefore, this loop checks all blocks which overlap with the current
    //   one by periodicity.
    for (pluint iOverlap = 0; iOverlap < localInfo.getPeriodicOverlapWithRemoteData().size();
         ++iOverlap)
    {
        Overlap3D overlap = localInfo.getPeriodicOverlapWithRemoteData()[iOverlap].overlap;
        if (contained(iX, iY, iZ, overlap.getOriginalCoordinates())) {
            plint overlapId = overlap.getOverlapId();
            foundId.push_back(overlapId);
            SmartBulk3D bulk(sparseBlock, envelopeWidth, overlapId);
            foundX.push_back(bulk.toLocalX(iX - overlap.getShiftX()));
            foundY.push_back(bulk.toLocalY(iY - overlap.getShiftY()));
            foundZ.push_back(bulk.toLocalZ(iZ - overlap.getShiftZ()));
        }
    }
    return hasBulkCell;
}

plint MultiBlockManagement3D::getRefinementLevel() const
{
    return refinementLevel;
}

void MultiBlockManagement3D::setRefinementLevel(plint newLevel)
{
    refinementLevel = newLevel;
}

void MultiBlockManagement3D::changeEnvelopeWidth(plint newEnvelopeWidth)
{
    envelopeWidth = newEnvelopeWidth;
    localInfo = LocalMultiBlockInfo3D(sparseBlock, getThreadAttribution(), envelopeWidth);
}

bool MultiBlockManagement3D::equivalentTo(MultiBlockManagement3D const &rhs) const
{
    std::map<plint, Box3D> const &bulks = sparseBlock.getBulks();
    std::map<plint, Box3D>::const_iterator it = bulks.begin();
    bool equalThreadAttribution = true;
    for (; it != bulks.end(); ++it) {
        plint blockId = it->first;
        if (threadAttribution->getMpiProcess(blockId)
                != rhs.threadAttribution->getMpiProcess(blockId)
            || threadAttribution->getLocalThreadId(blockId)
                   != rhs.threadAttribution->getLocalThreadId(blockId))
        {
            equalThreadAttribution = false;
            break;
        }
    }
    return sparseBlock.equals(rhs.sparseBlock) && equalThreadAttribution
           && refinementLevel == rhs.refinementLevel;
}

MultiBlockManagement3D scale(MultiBlockManagement3D const &originalManagement, plint relativeLevel)
{
    return MultiBlockManagement3D(
        scale(originalManagement.getSparseBlockStructure(), relativeLevel),
        originalManagement.getThreadAttribution().clone(), originalManagement.getEnvelopeWidth(),
        originalManagement.getRefinementLevel() + relativeLevel);
}

MultiBlockManagement3D intersect(
    MultiBlockManagement3D const &originalManagement, Box3D subDomain, bool crop)
{
    return MultiBlockManagement3D(
        intersect(originalManagement.getSparseBlockStructure(), subDomain, crop),
        originalManagement.getThreadAttribution().clone(), originalManagement.getEnvelopeWidth(),
        originalManagement.getRefinementLevel());
}

MultiBlockManagement3D intersect(
    MultiBlockManagement3D const &originalManagement, Box3D subDomain, Box3D newBoundingBox)
{
    return MultiBlockManagement3D(
        intersect(originalManagement.getSparseBlockStructure(), subDomain, newBoundingBox),
        originalManagement.getThreadAttribution().clone(), originalManagement.getEnvelopeWidth(),
        originalManagement.getRefinementLevel());
}

MultiBlockManagement3D intersect(
    MultiBlockManagement3D const &management1, MultiBlockManagement3D const &management2, bool crop)
{
    return MultiBlockManagement3D(
        intersect(
            management1.getSparseBlockStructure(), management2.getSparseBlockStructure(), crop),
        management1.getThreadAttribution().clone(), management1.getEnvelopeWidth(),
        management1.getRefinementLevel());
}

// TODO: Suspicious that this addedUniqueBulk is not used
MultiBlockManagement3D extend(
    MultiBlockManagement3D const &management, Box3D addedBulk, Box3D addedUniqueBulk)
{
    std::vector<plint> newIds;
    SparseBlockStructure3D resultStructure =
        extend(management.getSparseBlockStructure(), addedBulk, addedBulk, newIds);
    std::vector<plint> mpiProcesses(newIds.size()), localThreads(newIds.size());
    for (pluint iNew = 0; iNew < newIds.size(); ++iNew) {
        // Default-attribute the newly created blocks to the main process.
        mpiProcesses[iNew] = global::mpi().bossId();
        localThreads[iNew] = 0;
    }
    return MultiBlockManagement3D(
        resultStructure,
        management.getThreadAttribution().extend(newIds, mpiProcesses, localThreads),
        management.getEnvelopeWidth(), management.getRefinementLevel());
}

MultiBlockManagement3D except(MultiBlockManagement3D const &management, Box3D exceptedBlock)
{
    std::map<plint, std::vector<plint> > remappedIds;
    SparseBlockStructure3D resultStructure =
        except(management.getSparseBlockStructure(), exceptedBlock, remappedIds);
    return MultiBlockManagement3D(
        resultStructure,
        management.getThreadAttribution().merge(management.getThreadAttribution(), remappedIds),
        management.getEnvelopeWidth(), management.getRefinementLevel());
}

MultiBlockManagement3D block_union(
    MultiBlockManagement3D const &management1, MultiBlockManagement3D const &management2)
{
    std::map<plint, std::vector<plint> > remappedIds;
    SparseBlockStructure3D resultStructure = block_union(
        management1.getSparseBlockStructure(), management2.getSparseBlockStructure(), remappedIds);
    return MultiBlockManagement3D(
        resultStructure,
        management1.getThreadAttribution().merge(management2.getThreadAttribution(), remappedIds),
        management1.getEnvelopeWidth(), management1.getRefinementLevel());
}

MultiBlockManagement3D align(
    MultiBlockManagement3D const &originalManagement,
    MultiBlockManagement3D const &partnerManagement)
{
    std::vector<plint> newIds;
    std::map<plint, std::vector<plint> > remappedFromPartner;
    SparseBlockStructure3D resultStructure = alignDistribution3D(
        originalManagement.getSparseBlockStructure(), partnerManagement.getSparseBlockStructure(),
        newIds, remappedFromPartner);

    // 1. Parallelize the left-over blocks (the ones which don't overlap
    //    with partnerManagement) evenly.
    ExplicitThreadAttribution attribution;
    plint numBlocks = (plint)newIds.size();
    plint numProcs = global::mpi().getSize();
    plint iBlock = 0;
    for (plint iProc = 0; iProc < numProcs; ++iProc) {
        plint localNumBlocks = numBlocks / numProcs;
        if (iProc < numBlocks % numProcs) {
            ++localNumBlocks;
        }
        for (plint iLocal = 0; iLocal < localNumBlocks; ++iLocal) {
            attribution.addBlock(newIds[iBlock], iProc);
            ++iBlock;
        }
    }
    // 2. Merge remapped ids into the thread attribution, and return a
    //    corresponding MultiBlockManagement3D object.
    return MultiBlockManagement3D(
        resultStructure,
        attribution.merge(partnerManagement.getThreadAttribution(), remappedFromPartner),
        originalManagement.getEnvelopeWidth(), originalManagement.getRefinementLevel());
}

MultiBlockManagement3D align(
    std::vector<Box3D> const &originalDomain, MultiBlockManagement3D const &alignWith,
    plint envelopeWidth, plint refinementLevel, bool crop)
{
    Box3D bbox(alignWith.getBoundingBox());
    if (crop && !originalDomain.empty()) {
        bbox = originalDomain[0];
        for (pluint i = 1; i < originalDomain.size(); ++i) {
            bbox = bound(bbox, originalDomain[i]);
        }
    }
    SparseBlockStructure3D originalSparseBlock(bbox);
    for (plint i = 0; i < (plint)originalDomain.size(); ++i) {
        originalSparseBlock.addBlock(originalDomain[i], i);
    }
    MultiBlockManagement3D originalManagement(
        originalSparseBlock, defaultMultiBlockPolicy3D().getThreadAttribution(), envelopeWidth,
        refinementLevel);
    return align(originalManagement, alignWith);
}

MultiBlockManagement3D reparallelize(
    MultiBlockManagement3D const &management, plint blockLx, plint blockLy, plint blockLz)
{
    SparseBlockStructure3D resultStructure =
        reparallelize(management.getSparseBlockStructure(), blockLx, blockLy, blockLz);
    plint numBlocks = resultStructure.nextIncrementalId();
    plint numProcs = global::mpi().getSize();
    // Create a thread attribution from scratch, by partitioning the
    //   available blocks equally.
    ExplicitThreadAttribution *threadAttribution = new ExplicitThreadAttribution;
    plint iBlock = 0;
    for (plint iProc = 0; iProc < numProcs; ++iProc) {
        plint localNumBlocks = numBlocks / numProcs;
        if (iProc < numBlocks % numProcs) {
            ++localNumBlocks;
        }
        for (plint iLocal = 0; iLocal < localNumBlocks; ++iLocal) {
            threadAttribution->addBlock(iBlock, iProc);
            ++iBlock;
        }
    }

    return MultiBlockManagement3D(
        resultStructure, threadAttribution, management.getEnvelopeWidth(),
        management.getRefinementLevel());
}

SmartBulk3D::SmartBulk3D(MultiBlockManagement3D const &management, plint blockId) :
    sparseBlock(management.getSparseBlockStructure()), envelopeWidth(management.getEnvelopeWidth())
{
    sparseBlock.getBulk(blockId, bulk);
}

SmartBulk3D::SmartBulk3D(
    SparseBlockStructure3D const &sparseBlock_, plint envelopeWidth_, plint blockId) :
    sparseBlock(sparseBlock_), envelopeWidth(envelopeWidth_)
{
    sparseBlock.getBulk(blockId, bulk);
}

SmartBulk3D::SmartBulk3D(
    SparseBlockStructure3D const &sparseBlock_, plint envelopeWidth_, Box3D const &bulk_) :
    sparseBlock(sparseBlock_), envelopeWidth(envelopeWidth_), bulk(bulk_)
{ }

Box3D SmartBulk3D::getBulk() const
{
    return bulk;
}

Box3D SmartBulk3D::computeEnvelope() const
{
    return bulk.enlarge(envelopeWidth);
}

Box3D SmartBulk3D::computeNonPeriodicEnvelope() const
{
    Box3D boundingBox = sparseBlock.getBoundingBox();
    return Box3D(
        std::max(bulk.x0 - envelopeWidth, boundingBox.x0),
        std::min(bulk.x1 + envelopeWidth, boundingBox.x1),

        std::max(bulk.y0 - envelopeWidth, boundingBox.y0),
        std::min(bulk.y1 + envelopeWidth, boundingBox.y1),

        std::max(bulk.z0 - envelopeWidth, boundingBox.z0),
        std::min(bulk.z1 + envelopeWidth, boundingBox.z1));
}

Box3D SmartBulk3D::toLocal(Box3D const &coord) const
{
    return Box3D(
        coord.shift(-bulk.x0 + envelopeWidth, -bulk.y0 + envelopeWidth, -bulk.z0 + envelopeWidth));
}

plint SmartBulk3D::toLocalX(plint iX) const
{
    return iX - bulk.x0 + envelopeWidth;
}

plint SmartBulk3D::toLocalY(plint iY) const
{
    return iY - bulk.y0 + envelopeWidth;
}

plint SmartBulk3D::toLocalZ(plint iZ) const
{
    return iZ - bulk.z0 + envelopeWidth;
}

}  // namespace plb
