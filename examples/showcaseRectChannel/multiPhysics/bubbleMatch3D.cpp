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

#include "multiPhysics/bubbleMatch3D.h"

#include <limits>

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataField3D.hh"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/dataProcessingFunctional3D.hh"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "dataProcessors/dataInitializerFunctional3D.h"
#include "dataProcessors/dataInitializerFunctional3D.hh"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "dataProcessors/dataInitializerWrapper3D.hh"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/multiBlockGenerator3D.hh"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/multiDataField3D.hh"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.hh"
#include "multiBlock/nonLocalTransfer3D.h"
#include "multiBlock/nonLocalTransfer3D.hh"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/serialMultiDataField3D.h"
#include "multiBlock/serialMultiDataField3D.hh"
#include "multiPhysics/freeSurfaceUtil3D.h"
#include "offLattice/makeSparse3D.h"
#include "parallelism/mpiManager.h"
#include "parallelism/parallelMultiDataField3D.h"
#include "parallelism/parallelMultiDataField3D.hh"

namespace plb {

/* ************** class BubbleMatch3D ********************************** */

BubbleMatch3D::BubbleMatch3D(MultiBlock3D &templ) :
    bubbleContainer(createContainerBlock(templ, new BubbleCounterData3D(maxNumBubbles))),
    bubbleAnalysisContainer(createContainerBlock(templ, new BubbleAnalysisData3D())),
    bubbleRemapContainer(createContainerBlock(templ, new BubbleRemapData3D(maxNumBubbles))),
    mpiData(*bubbleContainer),
    tagMatrix(new MultiScalarField3D<plint>(*bubbleContainer))
{
    setToConstant(*tagMatrix, tagMatrix->getBoundingBox(), (plint)-1);
}

BubbleMatch3D::~BubbleMatch3D()
{
    delete bubbleContainer;
    delete bubbleAnalysisContainer;
    delete bubbleRemapContainer;
    delete tagMatrix;
}

pluint BubbleMatch3D::countAndTagBubbles()
{
    std::vector<MultiBlock3D *> args;
    args.push_back(tagMatrix);
    args.push_back(bubbleRemapContainer);
    applyProcessingFunctional(
        new CollectBubbleTags3D(), bubbleRemapContainer->getBoundingBox(), args);
    plint numBubbles = globalBubbleIds();
    applyProcessingFunctional(new ApplyTagRemap3D(), bubbleRemapContainer->getBoundingBox(), args);
    return numBubbles;
}

void BubbleMatch3D::computeBubbleData(pluint numBubbles)
{
    std::vector<double> bubbleCenterX(numBubbles), bubbleCenterY(numBubbles),
        bubbleCenterZ(numBubbles);
    bubbleVolume.resize(numBubbles);
    bubbleCenter.resize(numBubbles);

    std::vector<plint> const &localIds = mpiData.getLocalIds();
    for (pluint i = 0; i < localIds.size(); ++i) {
        plint id = localIds[i];
        AtomicContainerBlock3D &atomicDataContainer = bubbleAnalysisContainer->getComponent(id);
        BubbleAnalysisData3D *pData =
            dynamic_cast<BubbleAnalysisData3D *>(atomicDataContainer.getData());
        PLB_ASSERT(pData);
        BubbleAnalysisData3D &data = *pData;

        std::vector<double> const &nextVolume = data.bubbleVolume;
        std::vector<Array<double, 3> > const &nextCenter = data.bubbleCenter;
        PLB_ASSERT(nextVolume.size() == numBubbles);
        PLB_ASSERT(nextCenter.size() == numBubbles);

        for (pluint i = 0; i < numBubbles; ++i) {
            bubbleVolume[i] += nextVolume[i];
            bubbleCenterX[i] += nextCenter[i][0];
            bubbleCenterY[i] += nextCenter[i][1];
            bubbleCenterZ[i] += nextCenter[i][2];
        }
    }

#ifdef PLB_MPI_PARALLEL
    global::mpi().allReduceVect(bubbleVolume, MPI_SUM);
    global::mpi().allReduceVect(bubbleCenterX, MPI_SUM);
    global::mpi().allReduceVect(bubbleCenterY, MPI_SUM);
    global::mpi().allReduceVect(bubbleCenterZ, MPI_SUM);
#endif

    static const double epsilon = std::numeric_limits<double>::epsilon() * 1.e4;
    for (pluint i = 0; i < numBubbles; ++i) {
        bubbleCenter[i] = Array<double, 3>(bubbleCenterX[i], bubbleCenterY[i], bubbleCenterZ[i]);
        double volume = bubbleVolume[i];
        if (volume > epsilon) {
            bubbleCenter[i] /= volume;
        }
    }
}

plint BubbleMatch3D::globalBubbleIds()
{
    plint localNumUniqueBubbles = 0;
    std::vector<plint> const &localIds = mpiData.getLocalIds();
    for (pluint i = 0; i < localIds.size(); ++i) {
        plint id = localIds[i];
        AtomicContainerBlock3D &atomicDataContainer = bubbleRemapContainer->getComponent(id);
        BubbleRemapData3D *pData = dynamic_cast<BubbleRemapData3D *>(atomicDataContainer.getData());
        PLB_ASSERT(pData);
        BubbleRemapData3D &data = *pData;
        localNumUniqueBubbles += data.getUniqueTags().size();
    }

    std::vector<plint> allNumBubbles(global::mpi().getSize());
    allNumBubbles[global::mpi().getRank()] = localNumUniqueBubbles;
#ifdef PLB_MPI_PARALLEL
    global::mpi().allReduceVect(allNumBubbles, MPI_SUM);
#endif

    std::vector<plint> cumNumBubbles(global::mpi().getSize());
    PLB_ASSERT(cumNumBubbles.size() > 0);
    cumNumBubbles[0] = allNumBubbles[0];
    for (pluint i = 1; i < cumNumBubbles.size(); ++i) {
        cumNumBubbles[i] = allNumBubbles[i] + cumNumBubbles[i - 1];
    }
    plint totNumBubbles = cumNumBubbles.back();

    std::vector<plint> bubbleIds(totNumBubbles);
    pluint offset = global::mpi().getRank() == 0 ? 0 : cumNumBubbles[global::mpi().getRank() - 1];

    for (pluint i = 0; i < localIds.size(); ++i) {
        plint id = localIds[i];
        AtomicContainerBlock3D &atomicDataContainer = bubbleRemapContainer->getComponent(id);
        BubbleRemapData3D *pData = dynamic_cast<BubbleRemapData3D *>(atomicDataContainer.getData());
        PLB_ASSERT(pData);
        BubbleRemapData3D &data = *pData;
        std::vector<plint> uniqueTags = data.getUniqueTags();
        for (pluint i = 0; i < uniqueTags.size(); ++i, ++offset) {
            bubbleIds[offset] = uniqueTags[i];
        }
    }
#ifdef PLB_MPI_PARALLEL
    global::mpi().allReduceVect(bubbleIds, MPI_SUM);
#endif
    std::map<plint, plint> tagRemap;
    for (pluint i = 0; i < bubbleIds.size(); ++i) {
        tagRemap[bubbleIds[i]] = i;
    }

    for (pluint i = 0; i < localIds.size(); ++i) {
        plint id = localIds[i];
        AtomicContainerBlock3D &atomicDataContainer = bubbleRemapContainer->getComponent(id);
        BubbleRemapData3D *pData = dynamic_cast<BubbleRemapData3D *>(atomicDataContainer.getData());
        PLB_ASSERT(pData);
        BubbleRemapData3D &data = *pData;
        data.getTagRemap() = tagRemap;
    }

    return totNumBubbles;
}

void BubbleMatch3D::resetBubbleContainer()
{
    MultiBlockManagement3D const &management = bubbleContainer->getMultiBlockManagement();
    ThreadAttribution const &threadAttribution = management.getThreadAttribution();
    SparseBlockStructure3D const &sparseBlock = management.getSparseBlockStructure();
    std::map<plint, Box3D> const &domains = sparseBlock.getBulks();

    std::map<plint, Box3D>::const_iterator it = domains.begin();
    plint pos = 0;
    for (; it != domains.end(); ++it, ++pos) {
        plint id = it->first;
        if (threadAttribution.isLocal(id)) {
            AtomicContainerBlock3D &atomicDataContainer = bubbleContainer->getComponent(id);
            dynamic_cast<BubbleCounterData3D *>(atomicDataContainer.getData())->reset();
        }
    }
}

void BubbleMatch3D::bubbleBucketFill(MultiScalarField3D<int> &flag)
{
    setToConstant(*tagMatrix, tagMatrix->getBoundingBox(), (plint)-1);
    resetBubbleContainer();
    plint numIter = 2;
    while (numIter > 0) {
        CountBubbleIteration3D functional;
        std::vector<MultiBlock3D *> args;
        args.push_back(tagMatrix);
        args.push_back(&flag);
        args.push_back(bubbleContainer);
        applyProcessingFunctional(functional, bubbleContainer->getBoundingBox(), args);
        plint numConflicts = functional.getNumConflicts();
        if (numConflicts > 0) {
            numIter = 2;
        } else {
            --numIter;
        }
    }
}

/* ************** class BubbleMPIdata ********************************** */

BubbleMPIdata::BubbleMPIdata(MultiBlock3D &block)
{
    computeLocalIds(block);
}

std::vector<plint> const &BubbleMPIdata::getLocalIds() const
{
    return localIds;
}

void BubbleMPIdata::computeLocalIds(MultiBlock3D &block)
{
    MultiBlockManagement3D const &management = block.getMultiBlockManagement();
    ThreadAttribution const &threadAttribution = management.getThreadAttribution();
    SparseBlockStructure3D const &sparseBlock = management.getSparseBlockStructure();
    std::map<plint, Box3D> const &domains = sparseBlock.getBulks();

    std::map<plint, Box3D>::const_iterator it = domains.begin();
    plint pos = 0;
    for (; it != domains.end(); ++it, ++pos) {
        plint id = it->first;
        if (threadAttribution.isLocal(id)) {
            localIds.push_back(id);
        }
    }
}

/* ************** class BubbleCounterData3D ********************************** */

BubbleCounterData3D::BubbleCounterData3D(plint maxNumBubbles_) : maxNumBubbles(maxNumBubbles_) { }

BubbleCounterData3D *BubbleCounterData3D::clone() const
{
    return new BubbleCounterData3D(*this);
}

bool BubbleCounterData3D::convertCell(
    plint &tag0, plint tag1, plint tag2, plint tag3, plint tag4, plint tag5, plint tag6, plint tag7,
    plint tag8, plint tag9, plint tag10, plint tag11, plint tag12, plint tag13, plint tag1_,
    plint tag2_, plint tag3_, plint tag4_, plint tag5_, plint tag6_, plint tag7_, plint tag8_,
    plint tag9_, plint tag10_, plint tag11_, plint tag12_, plint tag13_)
{
    bool hasConflict = processNeighbor(tag0, tag1) || processNeighbor(tag0, tag2)
                       || processNeighbor(tag0, tag3) || processNeighbor(tag0, tag4)
                       || processNeighbor(tag0, tag5) || processNeighbor(tag0, tag6)
                       || processNeighbor(tag0, tag7) || processNeighbor(tag0, tag8)
                       || processNeighbor(tag0, tag9) || processNeighbor(tag0, tag10)
                       || processNeighbor(tag0, tag11) || processNeighbor(tag0, tag12)
                       || processNeighbor(tag0, tag13) || processNeighbor(tag0, tag1_)
                       || processNeighbor(tag0, tag2_) || processNeighbor(tag0, tag3_)
                       || processNeighbor(tag0, tag4_) || processNeighbor(tag0, tag5_)
                       || processNeighbor(tag0, tag6_) || processNeighbor(tag0, tag7_)
                       || processNeighbor(tag0, tag8_) || processNeighbor(tag0, tag9_)
                       || processNeighbor(tag0, tag10_) || processNeighbor(tag0, tag11_)
                       || processNeighbor(tag0, tag12_) || processNeighbor(tag0, tag13_);
    if (tag0 == -1) {
        tag0 = getNextTag();
    }
    return hasConflict;
}

bool BubbleCounterData3D::processNeighbor(plint &tag0, plint tag1)
{
    plint myTag = convertTag(tag0);
    plint otherTag = convertTag(tag1);
    tag0 = myTag;  // re-assign tag in case it got re-assigned in "convertTag()".
    bool hasConflict = false;
    if (otherTag == -1) {
        if (myTag == -1) {
            tag0 = getNextTag();
        }
    } else {
        if (myTag < otherTag) {
            tag0 = otherTag;
            registerConflict(myTag, otherTag);
            hasConflict = true;
        }
    }
    return hasConflict;
}

plint BubbleCounterData3D::getNextTag()
{
    plint nextTag = getUniqueID() * maxNumBubbles + nextCellId;
    ++nextCellId;
    return nextTag;
}

void BubbleCounterData3D::reset()
{
    nextCellId = 0;
    retagging.clear();
    // uniqueTags.clear();
    // tagRemap.clear();
}

plint BubbleCounterData3D::convertTag(plint tag) const
{
    if (tag == -1)
        return tag;
    std::map<plint, plint>::const_iterator it = retagging.find(tag);
    if (it == retagging.end()) {
        return tag;
    } else {
        return it->second;
    }
}

void BubbleCounterData3D::registerConflict(plint oldTag, plint newTag)
{
    retagging.insert(std::pair<plint, plint>(oldTag, newTag));
    std::map<plint, plint>::iterator it = retagging.begin();
    // If some of the map's items are still pointing to the old
    // tag, redirect them to the new one.
    for (; it != retagging.end(); ++it) {
        if (it->second == oldTag) {
            it->second = newTag;
        }
    }
}

/* *************** Class BubbleRemapData3D ******************************** */

BubbleRemapData3D *BubbleRemapData3D::clone() const
{
    return new BubbleRemapData3D(*this);
}

bool BubbleRemapData3D::isMyTag(plint tag)
{
    return tag / maxNumBubbles == getUniqueID();
}

/* *************** Class CountBubbleIteration3D ******************************** */

CountBubbleIteration3D::CountBubbleIteration3D() :
    numConflictsId(this->getStatistics().subscribeIntSum())
{ }

CountBubbleIteration3D *CountBubbleIteration3D::clone() const
{
    return new CountBubbleIteration3D(*this);
}

plint CountBubbleIteration3D::getNumConflicts() const
{
    return this->getStatistics().getIntSum(numConflictsId);
}

void CountBubbleIteration3D::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 3);
    ScalarField3D<plint> *pTagMatrix = dynamic_cast<ScalarField3D<plint> *>(atomicBlocks[0]);
    PLB_ASSERT(pTagMatrix);
    ScalarField3D<plint> &tagMatrix = *pTagMatrix;

    ScalarField3D<int> *pFlagMatrix = dynamic_cast<ScalarField3D<int> *>(atomicBlocks[1]);
    PLB_ASSERT(pFlagMatrix);
    ScalarField3D<int> &flagMatrix = *pFlagMatrix;

    AtomicContainerBlock3D *pDataBlock = dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[2]);
    PLB_ASSERT(pDataBlock);
    AtomicContainerBlock3D &dataBlock = *pDataBlock;
    BubbleCounterData3D *pData = dynamic_cast<BubbleCounterData3D *>(dataBlock.getData());
    PLB_ASSERT(pData);
    BubbleCounterData3D &data = *pData;

    Dot3D flagOffset = computeRelativeDisplacement(tagMatrix, flagMatrix);
    BlockStatistics &statistics = this->getStatistics();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                int currentFlag =
                    flagMatrix.get(iX + flagOffset.x, iY + flagOffset.y, iZ + flagOffset.z);
                if (currentFlag == freeSurfaceFlag::empty
                    || currentFlag == freeSurfaceFlag::interface) {
                    bool isConflicting = data.convertCell(
                        tagMatrix.get(iX, iY, iZ), tagMatrix.get(iX - 1, iY, iZ),
                        tagMatrix.get(iX, iY - 1, iZ), tagMatrix.get(iX, iY, iZ - 1),
                        tagMatrix.get(iX - 1, iY - 1, iZ), tagMatrix.get(iX - 1, iY + 1, iZ),
                        tagMatrix.get(iX - 1, iY, iZ - 1), tagMatrix.get(iX - 1, iY, iZ + 1),
                        tagMatrix.get(iX, iY - 1, iZ - 1), tagMatrix.get(iX, iY - 1, iZ + 1),
                        tagMatrix.get(iX - 1, iY - 1, iZ - 1),
                        tagMatrix.get(iX - 1, iY - 1, iZ + 1),
                        tagMatrix.get(iX - 1, iY + 1, iZ - 1),
                        tagMatrix.get(iX - 1, iY + 1, iZ + 1),

                        tagMatrix.get(iX + 1, iY, iZ), tagMatrix.get(iX, iY + 1, iZ),
                        tagMatrix.get(iX, iY, iZ + 1), tagMatrix.get(iX + 1, iY + 1, iZ),
                        tagMatrix.get(iX + 1, iY - 1, iZ), tagMatrix.get(iX + 1, iY, iZ + 1),
                        tagMatrix.get(iX + 1, iY, iZ - 1), tagMatrix.get(iX, iY + 1, iZ + 1),
                        tagMatrix.get(iX, iY + 1, iZ - 1), tagMatrix.get(iX + 1, iY + 1, iZ + 1),
                        tagMatrix.get(iX + 1, iY + 1, iZ - 1),
                        tagMatrix.get(iX + 1, iY - 1, iZ + 1),
                        tagMatrix.get(iX + 1, iY - 1, iZ - 1));
                    if (isConflicting) {
                        statistics.gatherIntSum(numConflictsId, 1);
                    }
                }
            }
        }
    }
}

/* *************** Class CollectBubbleTags3D ******************************** */

CollectBubbleTags3D *CollectBubbleTags3D::clone() const
{
    return new CollectBubbleTags3D(*this);
}

void CollectBubbleTags3D::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 2);
    ScalarField3D<plint> *pTagMatrix = dynamic_cast<ScalarField3D<plint> *>(atomicBlocks[0]);
    PLB_ASSERT(pTagMatrix);
    ScalarField3D<plint> &tagMatrix = *pTagMatrix;

    AtomicContainerBlock3D *pDataBlock = dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[1]);
    PLB_ASSERT(pDataBlock);
    AtomicContainerBlock3D &dataBlock = *pDataBlock;
    BubbleRemapData3D *pData = dynamic_cast<BubbleRemapData3D *>(dataBlock.getData());
    PLB_ASSERT(pData);
    BubbleRemapData3D &data = *pData;

    std::set<plint> uniqueTags;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint tag = tagMatrix.get(iX, iY, iZ);
                if (tag >= 0) {
                    uniqueTags.insert(tag);
                }
            }
        }
    }

    data.getUniqueTags().clear();
    std::set<plint>::const_iterator it = uniqueTags.begin();
    for (; it != uniqueTags.end(); ++it) {
        plint tag = *it;
        if (data.isMyTag(tag)) {
            data.getUniqueTags().push_back(tag);
        }
    }
}

/* *************** Class ApplyTagRemap3D ******************************** */

ApplyTagRemap3D *ApplyTagRemap3D::clone() const
{
    return new ApplyTagRemap3D(*this);
}

void ApplyTagRemap3D::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 2);
    ScalarField3D<plint> *pTagMatrix = dynamic_cast<ScalarField3D<plint> *>(atomicBlocks[0]);
    PLB_ASSERT(pTagMatrix);
    ScalarField3D<plint> &tagMatrix = *pTagMatrix;

    AtomicContainerBlock3D *pDataBlock = dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[1]);
    PLB_ASSERT(pDataBlock);
    AtomicContainerBlock3D &dataBlock = *pDataBlock;
    BubbleRemapData3D *pData = dynamic_cast<BubbleRemapData3D *>(dataBlock.getData());
    PLB_ASSERT(pData);
    BubbleRemapData3D &data = *pData;

    std::map<plint, plint> const &tagRemap = data.getTagRemap();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint tag = tagMatrix.get(iX, iY, iZ);
                if (tag >= 0) {
                    std::map<plint, plint>::const_iterator it = tagRemap.find(tag);
                    PLB_ASSERT(it != tagRemap.end());
                    tagMatrix.get(iX, iY, iZ) = it->second;
                }
            }
        }
    }
}

}  // namespace plb
