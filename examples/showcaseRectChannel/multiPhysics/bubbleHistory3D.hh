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

#ifndef BUBBLE_HISTORY_3D_HH
#define BUBBLE_HISTORY_3D_HH

#include <limits>

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "multiPhysics/bubbleHistory3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"
#include "offLattice/makeSparse3D.h"
#include "parallelism/mpiManager.h"

namespace plb {

/* ************** class BubbleHistory3D ********************************** */

template <typename T>
BubbleHistory3D<T>::BubbleHistory3D(MultiBlock3D &templ) :
    bubbleAnalysisContainer(createContainerBlock(templ, new BubbleAnalysisData3D())),
    bubbleCorrelationContainer(createContainerBlock(templ, new BubbleCorrelationData3D())),
    bubbleRemapContainer(createContainerBlock(templ, new BubbleRemapData3D())),
    mpiData(templ),
    oldTagMatrix(new MultiScalarField3D<plint>(templ)),
    nextBubbleID(0)
{
    setToConstant(*oldTagMatrix, oldTagMatrix->getBoundingBox(), (plint)-1);
}

template <typename T>
BubbleHistory3D<T>::~BubbleHistory3D()
{
    delete bubbleAnalysisContainer;
    delete bubbleCorrelationContainer;
    delete bubbleRemapContainer;
    delete oldTagMatrix;
}

template <typename T>
void BubbleHistory3D<T>::transition(
    BubbleMatch3D &bubbleMatch, plint iterationStep, T newBubbleVolumeCorrection,
    bool entrapBubbles, plint numRememberedVolumes)
{
    matchAndRemapBubbles(
        bubbleMatch, *oldTagMatrix, *bubbleMatch.getTagMatrix(), newBubbleVolumeCorrection,
        iterationStep, entrapBubbles, numRememberedVolumes);

    // Get the tag-matrix from bubbleMatch, and store it as internal oldTagMatrix for next
    // iteration.
    MultiScalarField3D<plint> *newTagMatrix = bubbleMatch.getTagMatrix();
    bubbleMatch.setTagMatrix(oldTagMatrix);
    oldTagMatrix = newTagMatrix;
}

template <typename T>
void BubbleHistory3D<T>::updateBubblePressure(
    MultiScalarField3D<T> &outsideDensity, T rhoEmpty, T alpha, T beta, T gamma)
{
    applyProcessingFunctional(
        new UpdateBubblePressure3D<T>(bubbles, rhoEmpty, alpha, beta, gamma),
        oldTagMatrix->getBoundingBox(), *oldTagMatrix, outsideDensity);
}

template <typename T>
void BubbleHistory3D<T>::freeze()
{
    typename std::map<plint, BubbleInfo3D>::iterator it = bubbles.begin();
    for (; it != bubbles.end(); ++it) {
        plint bubbleID = it->first;
        it->second.freeze();
        PLB_ASSERT((plint)fullBubbleRecord.size() > bubbleID);
        fullBubbleRecord[bubbleID].frozen = true;
    }
}

template <typename T>
void BubbleHistory3D<T>::freezeLargestBubble()
{
    typename std::map<plint, BubbleInfo3D>::iterator it;

    T biggestVolume = -1.0;
    plint largestBubbleID = -1;
    for (it = bubbles.begin(); it != bubbles.end(); ++it) {
        plint bubbleID = it->first;
        T volume = it->second.getVolume();
        if (volume > biggestVolume) {
            biggestVolume = volume;
            largestBubbleID = bubbleID;
        }
    }

    if (largestBubbleID != -1) {
        bubbles[largestBubbleID].freeze();
        PLB_ASSERT((plint)fullBubbleRecord.size() > largestBubbleID);
        fullBubbleRecord[largestBubbleID].frozen = true;
    }
}

template <typename T>
void BubbleHistory3D<T>::timeHistoryLog(std::string fName)
{
    plb_ofstream ofile(fName.c_str());
    std::map<plint, std::pair<std::vector<plint>, std::vector<plint> > >::const_iterator it =
        timeHistory.begin();
    for (; it != timeHistory.end(); ++it) {
        plint timeStep = it->first;
        std::pair<std::vector<plint>, std::vector<plint> > logEntry = it->second;
        ofile << "Iteration " << std::setw(8) << timeStep << ".  ";
        if (!logEntry.first.empty()) {
            ofile << "Created bubble(s) with ID";
            for (pluint i = 0; i < logEntry.first.size(); ++i) {
                ofile << " " << logEntry.first[i];
            }
            ofile << ".";
        }
        if (!logEntry.second.empty()) {
            ofile << " Removed bubble(s) with ID";
            for (pluint i = 0; i < logEntry.second.size(); ++i) {
                ofile << " " << logEntry.second[i];
            }
            ofile << ".";
        }
        ofile << std::endl;
    }
}

template <typename T>
void BubbleHistory3D<T>::fullBubbleLog(std::string fName)
{
    plb_ofstream ofile(fName.c_str());
    for (pluint i = 0; i < fullBubbleRecord.size(); ++i) {
        ofile << fullBubbleRecord[i].description(i);
        ofile << std::endl;
    }
}

template <typename T>
void BubbleHistory3D<T>::matchAndRemapBubbles(
    BubbleMatch3D &bubbleMatch, MultiScalarField3D<plint> &tagMatrix1,
    MultiScalarField3D<plint> &tagMatrix2, T newBubbleVolumeCorrection, plint iterationStep,
    bool entrapBubbles, plint numRememberedVolumes)
{
    pluint numBubbles = bubbleMatch.numBubbles();
    // The following data processor creates a correspondance between new bubbles and all overlapping
    // old bubbles. This is the basis of the upcoming contrived algorithm.
    std::vector<std::vector<plint> > newToAllOldMap;
    correlateBubbleIds(tagMatrix1, tagMatrix2, newToAllOldMap, numBubbles);

    // First, the map is inverted, to get a correspondence between old bubbles and all overlapping
    // new bubbles.
    std::map<plint, std::vector<plint> > oldToAllNewMap;
    typename std::map<plint, BubbleInfo3D>::const_iterator itOld = bubbles.begin();
    for (; itOld != bubbles.end(); ++itOld) {
        plint oldID = itOld->first;
        oldToAllNewMap.insert(std::pair<plint, std::vector<plint> >(oldID, std::vector<plint>()));
    }
    for (plint newID = 0; newID < (plint)newToAllOldMap.size(); ++newID) {
        for (pluint i = 0; i < newToAllOldMap[newID].size(); ++i) {
            plint oldID = newToAllOldMap[newID][i];
            typename std::map<plint, std::vector<plint> >::iterator it = oldToAllNewMap.find(oldID);
            PLB_ASSERT(it != oldToAllNewMap.end());
            it->second.push_back(newID);
        }
    }

    // Then, a transition pattern between groups of old and new bubbles is determined.
    std::vector<BubbleTransition3D> bubbleTransitions;
    computeBubbleTransitions(newToAllOldMap, oldToAllNewMap, bubbleTransitions);

    std::vector<plint> bubbleFinalRemap(newToAllOldMap.size());
    plint maxInt = std::numeric_limits<plint>::max();
    for (pluint i = 0; i < bubbleFinalRemap.size(); ++i) {
        plint minTag = maxInt;
        for (pluint j = 0; j < newToAllOldMap[i].size(); ++j) {
            if (newToAllOldMap[i][j] < minTag) {
                minTag = newToAllOldMap[i][j];
            }
        }
        bubbleFinalRemap[i] = minTag;
    }

    updateBubbleInformation(
        bubbleTransitions, bubbleMatch.getBubbleVolume(), newBubbleVolumeCorrection, iterationStep,
        entrapBubbles, numRememberedVolumes);

    std::vector<MultiBlock3D *> args2;
    args2.push_back(&tagMatrix2);
    args2.push_back(bubbleRemapContainer);
    applyProcessingFunctional(new ApplyTagRemap3D(), bubbleRemapContainer->getBoundingBox(), args2);
}

template <typename T>
void BubbleHistory3D<T>::updateBubbleInformation(
    std::vector<BubbleTransition3D> &bubbleTransitions, std::vector<double> const &bubbleVolume,
    T newBubbleVolumeCorrection, plint iterationStep, bool entrapBubbles,
    plint numRememberedVolumes)
{
    std::map<plint, BubbleInfo3D> newBubbles;
    std::map<plint, plint> newToFinal;
    for (pluint i = 0; i < bubbleTransitions.size(); ++i) {
        std::set<plint> &oldIDs = bubbleTransitions[i].oldIDs;
        std::set<plint> &newIDs = bubbleTransitions[i].newIDs;
        computeNewBubbles(
            oldIDs, newIDs, bubbleVolume, newBubbleVolumeCorrection, newBubbles, newToFinal,
            entrapBubbles, numRememberedVolumes);
        updateBubbleLog(bubbleTransitions[i], bubbleVolume, iterationStep, newBubbles, newToFinal);
    }

    // Prepare the data for remapping the bubbles.
    std::vector<plint> const &localIds = mpiData.getLocalIds();
    for (pluint i = 0; i < localIds.size(); ++i) {
        plint id = localIds[i];
        AtomicContainerBlock3D &atomicDataContainer = bubbleRemapContainer->getComponent(id);
        BubbleRemapData3D *pData = dynamic_cast<BubbleRemapData3D *>(atomicDataContainer.getData());
        PLB_ASSERT(pData);
        BubbleRemapData3D &data = *pData;
        data.getTagRemap() = newToFinal;
    }
    bubbles.swap(newBubbles);
}

template <typename T>
void BubbleHistory3D<T>::computeNewBubbles(
    std::set<plint> &oldIDs, std::set<plint> &newIDs, std::vector<double> const &bubbleVolume,
    T newBubbleVolumeCorrection, std::map<plint, BubbleInfo3D> &newBubbles,
    std::map<plint, plint> &newToFinal, bool entrapBubbles, plint numRememberedVolumes)
{
    // In this function there is code which implements two additional ideas:
    // 1.] Bubbles with a volume smaller than the cutoff-volume get automatically frozen.
    //     The cutoff value is given by the getCutoffVolume() function (this feature
    //     can be deactivated if the cutoff value is set to 0).
    // 2.] In the case of "bubble entrapment", if an original bubble is frozen and it splits
    //     into new bubbles, then only the largest of the new bubbles remains frozen.

    // Cutoff volume for freezing very small bubbles.
    static const T cutoffVolume = getCutoffVolume();

    // A bubble is created from nothing.
    if (oldIDs.empty()) {
        // Create a new bubble.
        PLB_ASSERT(newIDs.size() == 1);
        plint newID = *newIDs.begin();
        newBubbles[nextBubbleID] =
            BubbleInfo3D(bubbleVolume[newID] * newBubbleVolumeCorrection, numRememberedVolumes);
        // If a newly created bubble is smaller than the cutoff value, it gets automatically frozen.
        if (newBubbles[nextBubbleID].getVolume() < cutoffVolume) {
            newBubbles[nextBubbleID].freeze();
        }
        newToFinal[newID] = nextBubbleID;
        ++nextBubbleID;
    }
    // A bubble vanishes into nothing.
    else if (newIDs.empty())
    {
        PLB_ASSERT(oldIDs.size() == 1);
    }
    // All other cases: normal transition, splitting, merging.
    else
    {
        PLB_ASSERT(oldIDs.size() >= 1 && newIDs.size() >= 1);

        // First check the old bubbles: their total initial volume will be distributed
        // as initial volume of the new bubbles. If one of the old bubbles is frozen,
        // then so are the new bubbles.
        // And while we're at it, also close the history record for these old bubbles.
        T totalReferenceVolume = T();
        bool isFrozen = false;
        std::set<plint>::const_iterator itOld = oldIDs.begin();
        for (; itOld != oldIDs.end(); ++itOld) {
            plint oldID = *itOld;
            typename std::map<plint, BubbleInfo3D>::const_iterator it = bubbles.find(oldID);
            PLB_ASSERT(it != bubbles.end());
            totalReferenceVolume += it->second.getReferenceVolume();
            isFrozen = isFrozen || it->second.isFrozen();
        }

        // Then, compute the total volume of the new bubbles.
        // In the case of bubble entrapment, determine which
        // one is largest. Because, if the original bubble was
        // frozen, only the largest new one remains frozen.
        T newTotalVolume = T();
        std::set<plint>::const_iterator itNew = newIDs.begin();
        T largestVolume = (T)-1.0;
        std::set<plint>::const_iterator itLargest = newIDs.end();
        for (; itNew != newIDs.end(); ++itNew) {
            plint newID = *itNew;
            PLB_ASSERT(newID <= (plint)bubbleVolume.size());
            if (entrapBubbles && bubbleVolume[newID] > largestVolume) {
                largestVolume = bubbleVolume[newID];
                itLargest = itNew;
            }
            newTotalVolume += bubbleVolume[newID];
        }

        // Finally, put into place all information of the new bubbles.
        static const T epsilon = std::numeric_limits<T>::epsilon() * 1.e4;
        itNew = newIDs.begin();
        for (; itNew != newIDs.end(); ++itNew) {
            plint newID = *itNew;
            PLB_ASSERT(newID <= (plint)bubbleVolume.size());
            T newInitialVolume = bubbleVolume[newID];
            if (std::fabs(newTotalVolume) > epsilon) {
                newInitialVolume *= totalReferenceVolume / newTotalVolume;
            }

            plint finalID = newID;
            // In case of normal transition (no splitting nor merging), the
            // ID is inherited from the previous bubble.
            if (newIDs.size() == 1 && oldIDs.size() == 1) {
                finalID = *oldIDs.begin();
            }
            // Otherwise, a new ID is generated.
            else
            {
                finalID = nextBubbleID;
                ++nextBubbleID;
            }

            BubbleInfo3D nextBubble(newInitialVolume, numRememberedVolumes);
            nextBubble.setVolume(bubbleVolume[newID]);

            // If a bubble created through splitting/merging is smaller than
            // the cutoff value, it gets automatically frozen.
            if (nextBubble.getVolume() < cutoffVolume) {
                nextBubble.freeze();
            }

            if (!entrapBubbles) {
                if (isFrozen) {
                    nextBubble.freeze();
                }
            } else {
                // If original bubble was frozen, freeze the largest resulting bubble.
                if (isFrozen && itNew == itLargest) {
                    nextBubble.freeze();
                }
            }

            newToFinal[newID] = finalID;
            newBubbles[finalID] = nextBubble;
        }
    }
}

template <typename T>
void BubbleHistory3D<T>::updateBubbleLog(
    BubbleTransition3D &bubbleTransition, std::vector<double> const &bubbleVolume,
    plint iterationStep, std::map<plint, BubbleInfo3D> &newBubbles,
    std::map<plint, plint> &newToFinal)
{
    std::set<plint> &oldIDs = bubbleTransition.oldIDs;
    std::set<plint> &newIDs = bubbleTransition.newIDs;
    // Update the bubble-transition to point to the final ID (instead of
    // newID), wo it is properly written into the Log.
    std::set<plint> save_newIDs;
    for (std::set<plint>::iterator it = newIDs.begin(); it != newIDs.end(); ++it) {
        save_newIDs.insert(newToFinal[*it]);
    }
    save_newIDs.swap(newIDs);

    // For purposes of history, keep a full record of the bubble (start off the record).
    // Furthermore, keep a time-record of the fact that a bubble was created.
    //
    // A bubble is created from nothing.
    if (oldIDs.empty()) {
        PLB_ASSERT(newIDs.size() == 1);
        plint newID = *save_newIDs.begin();
        FullBubbleRecord nextBubbleRecord(bubbleVolume[newID], iterationStep);
        nextBubbleRecord.beginTransition = bubbleTransition;
        fullBubbleRecord.push_back(nextBubbleRecord);
        timeHistory[iterationStep].first.push_back(newToFinal[newID]);
    }
    // A bubble vanishes into nothing.
    else if (newIDs.empty())
    {
        PLB_ASSERT(oldIDs.size() == 1);
        plint vanishedID = *oldIDs.begin();
        PLB_ASSERT((plint)fullBubbleRecord.size() > vanishedID);
        fullBubbleRecord[vanishedID].endIteration = iterationStep;
        fullBubbleRecord[vanishedID].finalVolume = 0.;
        fullBubbleRecord[vanishedID].endTransition = bubbleTransition;
        timeHistory[iterationStep].second.push_back(vanishedID);
    }
    // Splitting or merging. Exclude the case of normal transitions.
    else if (!(newIDs.size() == 1 && oldIDs.size() == 1))
    {
        std::set<plint>::const_iterator itOld = oldIDs.begin();
        // Finish the record for all old bubbles that now vanish.
        for (; itOld != oldIDs.end(); ++itOld) {
            plint oldID = *itOld;
            typename std::map<plint, BubbleInfo3D>::const_iterator it = bubbles.find(oldID);
            PLB_ASSERT(it != bubbles.end());
            PLB_ASSERT((plint)fullBubbleRecord.size() > oldID);
            fullBubbleRecord[oldID].endIteration = iterationStep;
            fullBubbleRecord[oldID].finalVolume = it->second.getVolume();
            fullBubbleRecord[oldID].endTransition = bubbleTransition;
            timeHistory[iterationStep].second.push_back(oldID);
        }
        // Introduce a record for all newly created bubbles.
        std::set<plint>::const_iterator itNew = save_newIDs.begin();
        for (; itNew != save_newIDs.end(); ++itNew) {
            plint newID = *itNew;
            FullBubbleRecord nextBubbleRecord(bubbleVolume[newID], iterationStep);
            nextBubbleRecord.beginTransition = bubbleTransition;
            nextBubbleRecord.frozen = newBubbles[newToFinal[newID]].isFrozen();
            fullBubbleRecord.push_back(nextBubbleRecord);
            timeHistory[iterationStep].first.push_back(newToFinal[newID]);
        }
    }
}

template <typename T>
void BubbleHistory3D<T>::correlateBubbleIds(
    MultiScalarField3D<plint> &tagMatrix1, MultiScalarField3D<plint> &tagMatrix2,
    std::vector<std::vector<plint> > &newToAllOldMap, pluint numBubbles)
{
    std::vector<plint> const &localIds = mpiData.getLocalIds();
    plint maxInt = std::numeric_limits<plint>::max();

    // Reset the vector for the return value.
    newToAllOldMap = std::vector<std::vector<plint> >(numBubbles);

    // Initialization: prepare the new-to-old maps in the containers of the
    // local MPI threads.
    for (pluint i = 0; i < localIds.size(); ++i) {
        plint id = localIds[i];
        AtomicContainerBlock3D &atomicDataContainer = bubbleCorrelationContainer->getComponent(id);
        BubbleCorrelationData3D *pData =
            dynamic_cast<BubbleCorrelationData3D *>(atomicDataContainer.getData());
        PLB_ASSERT(pData);
        BubbleCorrelationData3D &data = *pData;

        data.newToOldMap0.resize(numBubbles);
        data.newToOldMap1.resize(numBubbles);
        // newToOldMap0 holds the "previous iteration": initialize it to -1
        // to make sure it is smaller than any actually occurring tag.
        std::fill(data.newToOldMap0.begin(), data.newToOldMap0.end(), -1);
        // newToOldMap1 will take the tags from the first iteration. Initialize
        // it to maxInt in view of the coming MIN reduction.
        std::fill(data.newToOldMap1.begin(), data.newToOldMap1.end(), maxInt);
    }

    // The number of iterations in the coming loop will be N+1, where N is the maximum
    // number of bubbles from the physical time step t-1 that coalesced to form a bubble
    // at time t. Very often, no coalescence occurs and N is just 1.
    bool hasConverged = false;
    while (!hasConverged) {
        std::vector<MultiBlock3D *> args;
        args.push_back(&tagMatrix1);
        args.push_back(&tagMatrix2);
        args.push_back(bubbleCorrelationContainer);
        // Correlate the new bubbles with the next group of old bubbles.
        applyProcessingFunctional(
            new CorrelateBubbleIds3D<T>, bubbleCorrelationContainer->getBoundingBox(), args);

        std::vector<plint> newToOldMap(numBubbles, maxInt);

        // Apply the MIN reduction to all atomic-blocks that are on the current MPI thread.
        for (pluint i = 0; i < localIds.size(); ++i) {
            plint id = localIds[i];
            AtomicContainerBlock3D &atomicDataContainer =
                bubbleCorrelationContainer->getComponent(id);
            BubbleCorrelationData3D *pData =
                dynamic_cast<BubbleCorrelationData3D *>(atomicDataContainer.getData());
            PLB_ASSERT(pData);
            BubbleCorrelationData3D &data = *pData;

            PLB_ASSERT(data.newToOldMap0.size() == numBubbles);
            PLB_ASSERT(data.newToOldMap1.size() == numBubbles);

            for (pluint i = 0; i < numBubbles; ++i) {
                newToOldMap[i] = std::min(newToOldMap[i], data.newToOldMap1[i]);
            }
        }

        // Apply the MIN reduction globally.
#ifdef PLB_MPI_PARALLEL
        global::mpi().allReduceVect(newToOldMap, MPI_MIN);
#endif

        // If new bubbles where found, add them to the return value. If not, abort the
        // iterations.
        bool newBubbleDiscovered = false;
        for (pluint i = 0; i < numBubbles; ++i) {
            if (newToOldMap[i] < maxInt) {
                newBubbleDiscovered = true;
                newToAllOldMap[i].push_back(newToOldMap[i]);
            }
        }
        hasConverged = !newBubbleDiscovered;

        // Assign the result of the global reduction to all local atomic-blocks, to be
        // used for the next iteration.
        for (pluint i = 0; i < localIds.size(); ++i) {
            plint id = localIds[i];
            AtomicContainerBlock3D &atomicDataContainer =
                bubbleCorrelationContainer->getComponent(id);
            BubbleCorrelationData3D *pData =
                dynamic_cast<BubbleCorrelationData3D *>(atomicDataContainer.getData());
            PLB_ASSERT(pData);
            BubbleCorrelationData3D &data = *pData;

            // Advance by one iteration: assign the globally reduced data to newToOldMap0
            // ("the previous iteration").
            data.newToOldMap0.assign(newToOldMap.begin(), newToOldMap.end());
            // Prepare newToOldMap1 for the next iteration, and especially for the coming
            // MIN reduction, by initializing to maxInt.
            std::fill(data.newToOldMap1.begin(), data.newToOldMap1.end(), maxInt);
        }
    }
}

/* *************** Class CorrelateBubbleIds3D ******************************** */

template <typename T>
CorrelateBubbleIds3D<T> *CorrelateBubbleIds3D<T>::clone() const
{
    return new CorrelateBubbleIds3D<T>(*this);
}

template <typename T>
void CorrelateBubbleIds3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 3);

    ScalarField3D<plint> *pTagMatrix1 = dynamic_cast<ScalarField3D<plint> *>(atomicBlocks[0]);
    PLB_ASSERT(pTagMatrix1);
    ScalarField3D<plint> &tagMatrix1 = *pTagMatrix1;

    ScalarField3D<plint> *pTagMatrix2 = dynamic_cast<ScalarField3D<plint> *>(atomicBlocks[1]);
    PLB_ASSERT(pTagMatrix2);
    ScalarField3D<plint> &tagMatrix2 = *pTagMatrix2;

    AtomicContainerBlock3D *pDataBlock = dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[2]);
    PLB_ASSERT(pDataBlock);
    AtomicContainerBlock3D &dataBlock = *pDataBlock;
    BubbleCorrelationData3D *pData = dynamic_cast<BubbleCorrelationData3D *>(dataBlock.getData());
    PLB_ASSERT(pData);
    BubbleCorrelationData3D &data = *pData;

    Dot3D tag2Ofs = computeRelativeDisplacement(tagMatrix1, tagMatrix2);

#ifdef PLB_DEBUG
    pluint numNewBubbles =
#endif
        data.newToOldMap0.size();
    PLB_ASSERT(data.newToOldMap1.size() == numNewBubbles);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint newBubbleTag = tagMatrix2.get(iX + tag2Ofs.x, iY + tag2Ofs.y, iZ + tag2Ofs.z);
                if (newBubbleTag >= 0) {
                    plint oldBubbleTag = tagMatrix1.get(iX, iY, iZ);
                    if (oldBubbleTag >= 0) {
                        PLB_ASSERT(newBubbleTag < (plint)numNewBubbles);
                        // At this step of the iteration (corresponding to a call to the present
                        // data processor), we're only targeting tags that are larger than the
                        // ones from the previous step.
                        if (oldBubbleTag > data.newToOldMap0[newBubbleTag]) {
                            // Searching the minimum among all tags that match the condition.
                            data.newToOldMap1[newBubbleTag] =
                                std::min(data.newToOldMap1[newBubbleTag], oldBubbleTag);
                        }
                    }
                }
            }
        }
    }
}

/* *************** Class UpdateBubblePressure3D ******************************** */

template <typename T>
UpdateBubblePressure3D<T>::UpdateBubblePressure3D(
    std::map<plint, BubbleInfo3D> const &bubbles_, T rho0_, T alpha_, T beta_, T gamma_) :
    bubbles(bubbles_), rho0(rho0_), alpha(alpha_), beta(beta_), gamma(gamma_)
{ }

template <typename T>
UpdateBubblePressure3D<T> *UpdateBubblePressure3D<T>::clone() const
{
    return new UpdateBubblePressure3D<T>(*this);
}

template <typename T>
void UpdateBubblePressure3D<T>::process(
    Box3D domain, ScalarField3D<plint> &tags, ScalarField3D<T> &density)
{
    // Cutoff density ratio for the bubble internal pressure.
    static const T cutoffDensityRatio = getCutoffDensityRatio();

    Dot3D ofs = computeRelativeDisplacement(tags, density);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint tag = tags.get(iX, iY, iZ);
                if (tag >= 0) {
                    typename std::map<plint, BubbleInfo3D>::const_iterator it = bubbles.find(tag);
                    PLB_ASSERT(it != bubbles.end());
                    BubbleInfo3D const &info = it->second;
                    T rho = rho0;
                    if (!info.isFrozen()) {
                        T volumeRatio = info.getVolumeRatio();
                        if (alpha < (T)0 || beta < (T)0) {
                            // Model: rho = rho0 * (V0/V)^gamma
                            rho = rho0 * std::pow(volumeRatio, gamma);
                        } else {
                            // Model: rho = rho0 * [1 + alpha * (1 - V/V0)^beta]
                            rho =
                                rho0 * (1.0 + alpha * std::pow((T)(1.0 - 1.0 / volumeRatio), beta));
                        }
                        if (rho > cutoffDensityRatio * rho0) {
                            rho = cutoffDensityRatio * rho0;
                        }
                    }
                    density.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z) = rho;
                } else {
                    density.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z) = rho0;
                }
            }
        }
    }
}

template <typename T>
void BubbleHistory3D<T>::computeBubbleTransitions(
    std::vector<std::vector<plint> > const &newToAllOldMap,
    std::map<plint, std::vector<plint> > const &oldToAllNewMap,
    std::vector<BubbleTransition3D> &bubbleTransitions)
{
    std::vector<bool> checked(newToAllOldMap.size());
    std::fill(checked.begin(), checked.end(), false);
    for (pluint newId = 0; newId < newToAllOldMap.size(); ++newId) {
        std::set<plint> newCandidates, oldCandidates;
        BubbleTransition3D transition;
        if (!checked[newId]) {
            newCandidates.insert(newId);
        }
        while (!(newCandidates.empty() && oldCandidates.empty())) {
            if (!newCandidates.empty()) {
                std::set<plint>::iterator it1 = newCandidates.begin();
                plint newCandidate = *it1;
                newCandidates.erase(it1);

                std::set<plint>::iterator it2 = transition.newIDs.find(newCandidate);
                // Treat this ID only if it is has not already been treated.
                if (it2 == transition.newIDs.end() && !checked[newCandidate]) {
                    transition.newIDs.insert(newCandidate);
                    checked[newCandidate] = true;

                    PLB_ASSERT(newCandidate < (plint)newToAllOldMap.size());
                    std::vector<plint> oldIds = newToAllOldMap[newCandidate];
                    oldCandidates.insert(oldIds.begin(), oldIds.end());
                }
            } else if (!oldCandidates.empty()) {
                std::set<plint>::iterator it1 = oldCandidates.begin();
                plint oldCandidate = *it1;
                oldCandidates.erase(it1);

                std::set<plint>::iterator it2 = transition.oldIDs.find(oldCandidate);
                // Treat this ID only if it is has not already been treated.
                if (it2 == transition.oldIDs.end()) {
                    transition.oldIDs.insert(oldCandidate);

                    std::map<plint, std::vector<plint> >::const_iterator it2 =
                        oldToAllNewMap.find(oldCandidate);
                    PLB_ASSERT(it2 != oldToAllNewMap.end());
                    std::vector<plint> newIds = it2->second;
                    newCandidates.insert(newIds.begin(), newIds.end());
                }
            }
        }
        if (!transition.empty()) {
            bubbleTransitions.push_back(transition);
        }
    }

    std::map<plint, std::vector<plint> >::const_iterator it = oldToAllNewMap.begin();
    for (; it != oldToAllNewMap.end(); ++it) {
        if (it->second.empty()) {
            BubbleTransition3D transition;
            transition.oldIDs.insert(it->first);
            bubbleTransitions.push_back(transition);
        }
    }
}

}  // namespace plb

#endif  // BUBBLE_HISTORY_3D_HH
