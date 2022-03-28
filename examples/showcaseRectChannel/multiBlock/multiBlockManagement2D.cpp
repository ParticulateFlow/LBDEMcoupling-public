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
 * Geometry specifications for 2D multiblocks -- implementation.
 */

#include "multiBlock/multiBlockManagement2D.h"

#include <algorithm>

#include "core/plbDebug.h"
#include "io/parallelIO.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiBlock/staticRepartitions2D.h"

namespace plb {

MultiBlockManagement2D::MultiBlockManagement2D(
    SparseBlockStructure2D const &sparseBlock_, ThreadAttribution *threadAttribution_,
    plint envelopeWidth_, plint refinementLevel_) :
    envelopeWidth(envelopeWidth_),
    sparseBlock(sparseBlock_),
    threadAttribution(threadAttribution_),
    localInfo(sparseBlock, getThreadAttribution(), envelopeWidth),
    refinementLevel(refinementLevel_)
{ }

MultiBlockManagement2D::MultiBlockManagement2D(MultiBlockManagement2D const &rhs) :
    envelopeWidth(rhs.envelopeWidth),
    sparseBlock(rhs.sparseBlock),
    threadAttribution(rhs.threadAttribution->clone()),
    localInfo(rhs.localInfo),
    refinementLevel(rhs.refinementLevel)
{ }

MultiBlockManagement2D &MultiBlockManagement2D::operator=(MultiBlockManagement2D const &rhs)
{
    MultiBlockManagement2D newBlockManagement(rhs);
    newBlockManagement.swap(*this);
    return *this;
}

void MultiBlockManagement2D::swap(MultiBlockManagement2D &rhs)
{
    std::swap(envelopeWidth, rhs.envelopeWidth);
    sparseBlock.swap(rhs.sparseBlock);
    std::swap(threadAttribution, rhs.threadAttribution);
    localInfo.swap(rhs.localInfo);
    std::swap(refinementLevel, rhs.refinementLevel);
}

MultiBlockManagement2D::~MultiBlockManagement2D()
{
    delete threadAttribution;
}

plint MultiBlockManagement2D::getEnvelopeWidth() const
{
    return envelopeWidth;
}

Box2D MultiBlockManagement2D::getBoundingBox() const
{
    return sparseBlock.getBoundingBox();
}

Box2D MultiBlockManagement2D::getBulk(plint blockId) const
{
    Box2D bulk;
#ifdef PLB_DEBUG
    bool ok =
#endif
        sparseBlock.getBulk(blockId, bulk);
    PLB_ASSERT(ok);
    return bulk;
}

Box2D MultiBlockManagement2D::getUniqueBulk(plint blockId) const
{
    Box2D uniqueBulk;
#ifdef PLB_DEBUG
    bool ok =
#endif
        sparseBlock.getUniqueBulk(blockId, uniqueBulk);
    PLB_ASSERT(ok);
    return uniqueBulk;
}

Box2D MultiBlockManagement2D::getEnvelope(plint blockId) const
{
    Box2D bulk;
#ifdef PLB_DEBUG
    bool ok =
#endif
        sparseBlock.getBulk(blockId, bulk);
    PLB_ASSERT(ok);
    return SmartBulk2D(sparseBlock, envelopeWidth, bulk).computeEnvelope();
}

SparseBlockStructure2D const &MultiBlockManagement2D::getSparseBlockStructure() const
{
    return sparseBlock;
}

LocalMultiBlockInfo2D const &MultiBlockManagement2D::getLocalInfo() const
{
    return localInfo;
}

ThreadAttribution const &MultiBlockManagement2D::getThreadAttribution() const
{
    return *threadAttribution;
}

bool MultiBlockManagement2D::findInLocalBulk(
    plint iX, plint iY, plint &foundId, plint &localX, plint &localY) const
{
    foundId = sparseBlock.locate(iX, iY);
    SmartBulk2D bulk(sparseBlock, envelopeWidth, foundId);
    localX = bulk.toLocalX(iX);
    localY = bulk.toLocalY(iY);
    return foundId >= 0;
}

bool MultiBlockManagement2D::findAllLocalRepresentations(
    plint iX, plint iY, std::vector<plint> &foundId, std::vector<plint> &foundX,
    std::vector<plint> &foundY) const
{
    bool hasBulkCell = false;
    // First, search in all blocks which are local to the current processor,
    // including in the envelopes.
    for (pluint iBlock = 0; iBlock < localInfo.getBlocks().size(); ++iBlock) {
        plint blockId = localInfo.getBlocks()[iBlock];
        SmartBulk2D bulk(sparseBlock, envelopeWidth, blockId);
        if (contained(iX, iY, bulk.computeEnvelope())) {
            if (contained(iX, iY, bulk.getBulk())) {
                hasBulkCell = true;
                foundId.insert(foundId.begin(), blockId);
                foundX.insert(foundX.begin(), bulk.toLocalX(iX));
                foundY.insert(foundY.begin(), bulk.toLocalY(iY));
            } else {
                foundId.push_back(blockId);
                foundX.push_back(bulk.toLocalX(iX));
                foundY.push_back(bulk.toLocalY(iY));
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
        Overlap2D overlap = localInfo.getPeriodicOverlapWithRemoteData()[iOverlap].overlap;
        if (contained(iX, iY, overlap.getOriginalCoordinates())) {
            plint overlapId = overlap.getOverlapId();
            foundId.push_back(overlapId);
            SmartBulk2D bulk(sparseBlock, envelopeWidth, overlapId);
            foundX.push_back(bulk.toLocalX(iX - overlap.getShiftX()));
            foundY.push_back(bulk.toLocalY(iY - overlap.getShiftY()));
        }
    }
    return hasBulkCell;
}

plint MultiBlockManagement2D::getRefinementLevel() const
{
    return refinementLevel;
}

void MultiBlockManagement2D::setRefinementLevel(plint newLevel)
{
    refinementLevel = newLevel;
}

void MultiBlockManagement2D::changeEnvelopeWidth(plint newEnvelopeWidth)
{
    envelopeWidth = newEnvelopeWidth;
    localInfo = LocalMultiBlockInfo2D(sparseBlock, getThreadAttribution(), envelopeWidth);
}

bool MultiBlockManagement2D::equivalentTo(MultiBlockManagement2D const &rhs) const
{
    std::map<plint, Box2D> const &bulks = sparseBlock.getBulks();
    std::map<plint, Box2D>::const_iterator it = bulks.begin();
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

MultiBlockManagement2D intersect(
    MultiBlockManagement2D const &originalManagement, Box2D subDomain, bool crop)
{
    return MultiBlockManagement2D(
        intersect(originalManagement.getSparseBlockStructure(), subDomain, crop),
        originalManagement.getThreadAttribution().clone(), originalManagement.getEnvelopeWidth(),
        originalManagement.getRefinementLevel());
}

MultiBlockManagement2D intersect(
    MultiBlockManagement2D const &originalManagement, Box2D subDomain, Box2D newBoundingBox)
{
    return MultiBlockManagement2D(
        intersect(originalManagement.getSparseBlockStructure(), subDomain, newBoundingBox),
        originalManagement.getThreadAttribution().clone(), originalManagement.getEnvelopeWidth(),
        originalManagement.getRefinementLevel());
}

MultiBlockManagement2D intersect(
    MultiBlockManagement2D const &management1, MultiBlockManagement2D const &management2, bool crop)
{
    return MultiBlockManagement2D(
        intersect(
            management1.getSparseBlockStructure(), management2.getSparseBlockStructure(), crop),
        management1.getThreadAttribution().clone(), management1.getEnvelopeWidth(),
        management1.getRefinementLevel());
}

// TODO: Suspicious that this unique bulk is not used
MultiBlockManagement2D extend(
    MultiBlockManagement2D const &management, Box2D addedBulk, Box2D addedUniqueBulk)
{
    std::vector<plint> newIds;
    SparseBlockStructure2D resultStructure =
        extend(management.getSparseBlockStructure(), addedBulk, addedBulk, newIds);
    std::vector<plint> mpiProcesses(newIds.size()), localThreads(newIds.size());
    for (pluint iNew = 0; iNew < newIds.size(); ++iNew) {
        // Default-attribute the newly created blocks to the main process.
        mpiProcesses[iNew] = global::mpi().bossId();
        localThreads[iNew] = 0;
    }
    return MultiBlockManagement2D(
        resultStructure,
        management.getThreadAttribution().extend(newIds, mpiProcesses, localThreads),
        management.getEnvelopeWidth(), management.getRefinementLevel());
}

MultiBlockManagement2D except(MultiBlockManagement2D const &management, Box2D exceptedBlock)
{
    std::map<plint, std::vector<plint> > remappedIds;
    SparseBlockStructure2D resultStructure =
        except(management.getSparseBlockStructure(), exceptedBlock, remappedIds);
    return MultiBlockManagement2D(
        resultStructure,
        management.getThreadAttribution().merge(management.getThreadAttribution(), remappedIds),
        management.getEnvelopeWidth(), management.getRefinementLevel());
}

MultiBlockManagement2D block_union(
    MultiBlockManagement2D const &management1, MultiBlockManagement2D const &management2)
{
    std::map<plint, std::vector<plint> > remappedIds;
    SparseBlockStructure2D resultStructure = block_union(
        management1.getSparseBlockStructure(), management2.getSparseBlockStructure(), remappedIds);
    return MultiBlockManagement2D(
        resultStructure,
        management1.getThreadAttribution().merge(management2.getThreadAttribution(), remappedIds),
        management1.getEnvelopeWidth(), management1.getRefinementLevel());
}

MultiBlockManagement2D align(
    MultiBlockManagement2D const &originalManagement,
    MultiBlockManagement2D const &partnerManagement)
{
    std::vector<plint> newIds;
    std::map<plint, std::vector<plint> > remappedFromPartner;
    SparseBlockStructure2D resultStructure = alignDistribution2D(
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
    //    corresponding MultiBlockManagement2D object.
    return MultiBlockManagement2D(
        resultStructure,
        attribution.merge(partnerManagement.getThreadAttribution(), remappedFromPartner),
        originalManagement.getEnvelopeWidth(), originalManagement.getRefinementLevel());
}

MultiBlockManagement2D align(
    std::vector<Box2D> const &originalDomain, MultiBlockManagement2D const &alignWith,
    plint envelopeWidth, plint refinementLevel, bool crop)
{
    Box2D bbox(alignWith.getBoundingBox());
    if (crop && !originalDomain.empty()) {
        bbox = originalDomain[0];
        for (pluint i = 1; i < originalDomain.size(); ++i) {
            bbox = bound(bbox, originalDomain[i]);
        }
    }
    SparseBlockStructure2D originalSparseBlock(bbox);
    for (plint i = 0; i < (plint)originalDomain.size(); ++i) {
        originalSparseBlock.addBlock(originalDomain[i], i);
    }
    MultiBlockManagement2D originalManagement(
        originalSparseBlock, defaultMultiBlockPolicy2D().getThreadAttribution(), envelopeWidth,
        refinementLevel);
    return align(originalManagement, alignWith);
}

MultiBlockManagement2D reparallelize(
    MultiBlockManagement2D const &management, plint blockLx, plint blockLy)
{
    SparseBlockStructure2D resultStructure =
        reparallelize(management.getSparseBlockStructure(), blockLx, blockLy);
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
    return MultiBlockManagement2D(
        resultStructure, threadAttribution, management.getEnvelopeWidth(),
        management.getRefinementLevel());
}

SmartBulk2D::SmartBulk2D(MultiBlockManagement2D const &management, plint blockId) :
    sparseBlock(management.getSparseBlockStructure()), envelopeWidth(management.getEnvelopeWidth())
{
    sparseBlock.getBulk(blockId, bulk);
}

SmartBulk2D::SmartBulk2D(
    SparseBlockStructure2D const &sparseBlock_, plint envelopeWidth_, plint blockId) :
    sparseBlock(sparseBlock_), envelopeWidth(envelopeWidth_)
{
    sparseBlock.getBulk(blockId, bulk);
}

SmartBulk2D::SmartBulk2D(
    SparseBlockStructure2D const &sparseBlock_, plint envelopeWidth_, Box2D const &bulk_) :
    sparseBlock(sparseBlock_), envelopeWidth(envelopeWidth_), bulk(bulk_)
{ }

Box2D SmartBulk2D::getBulk() const
{
    return bulk;
}

Box2D SmartBulk2D::computeEnvelope() const
{
    return bulk.enlarge(envelopeWidth);
}

Box2D SmartBulk2D::computeNonPeriodicEnvelope() const
{
    Box2D boundingBox = sparseBlock.getBoundingBox();
    return Box2D(
        std::max(bulk.x0 - envelopeWidth, boundingBox.x0),
        std::min(bulk.x1 + envelopeWidth, boundingBox.x1),

        std::max(bulk.y0 - envelopeWidth, boundingBox.y0),
        std::min(bulk.y1 + envelopeWidth, boundingBox.y1));
}

Box2D SmartBulk2D::toLocal(Box2D const &coord) const
{
    return Box2D(coord.shift(-bulk.x0 + envelopeWidth, -bulk.y0 + envelopeWidth));
}

plint SmartBulk2D::toLocalX(plint iX) const
{
    return iX - bulk.x0 + envelopeWidth;
}

plint SmartBulk2D::toLocalY(plint iY) const
{
    return iY - bulk.y0 + envelopeWidth;
}

}  // namespace plb
