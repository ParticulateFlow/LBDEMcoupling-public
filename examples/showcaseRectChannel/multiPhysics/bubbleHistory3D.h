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

#ifndef BUBBLE_HISTORY_3D_H
#define BUBBLE_HISTORY_3D_H

#include <sstream>
#include <string>

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "core/globalDefs.h"
#include "multiPhysics/bubbleMatch3D.h"
#include "multiPhysics/freeSurfaceModel3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"
#include "offLattice/makeSparse3D.h"
#include "parallelism/mpiManager.h"

namespace plb {

struct BubbleTransition3D;
class BubbleInfo3D;
struct FullBubbleRecord;

template <typename T>
class BubbleHistory3D {
public:
    BubbleHistory3D(MultiBlock3D &templ);
    ~BubbleHistory3D();
    void transition(
        BubbleMatch3D &bubbleMatch, plint iterationStep, T newBubbleVolumeCorrection = (T)1.,
        bool entrapBubbles = false, plint numRememberedVolumes = 1);
    // Based on the bubble information (original volume, current volume), update the pressure field
    // for all bubbles.
    void updateBubblePressure(
        MultiScalarField3D<T> &outsideDensity, T rhoEmpty, T alpha = (T)-1, T beta = (T)-1,
        T gamma = (T)1);
    void freeze();
    void freezeLargestBubble();
    void timeHistoryLog(std::string fName);
    void fullBubbleLog(std::string fName);
    std::map<plint, BubbleInfo3D> const &getBubbles() const
    {
        return bubbles;
    }

private:
    // Implement the time evolution of the bubbles:
    // - Reassign an ID to the new bubbles, compatible with the old ID.
    // - Create a data structure for the new bubbles.
    // - Print a message if bubbles are created/deleted.
    void matchAndRemapBubbles(
        BubbleMatch3D &bubbleMatch, MultiScalarField3D<plint> &tagMatrix1,
        MultiScalarField3D<plint> &tagMatrix2, T newBubbleVolumeCorrection, plint iterationStep,
        bool entrapBubbles, plint numRememberedVolumes);
    // Create the map "newToAllOldMap" which maps the temporary new IDs (continuously
    // numbered) to the "old" IDs (either existing ones or newly attributed).
    void correlateBubbleIds(
        MultiScalarField3D<plint> &tagMatrix1, MultiScalarField3D<plint> &tagMatrix2,
        std::vector<std::vector<plint> > &newToAllOldMap, pluint numBubbles);
    // Updates the (non-parallel) data structure which holds the overview of the currently
    // available bubbles.
    void updateBubbleInformation(
        std::vector<BubbleTransition3D> &bubbleTransitions, std::vector<double> const &bubbleVolume,
        T newBubbleVolumeCorrection, plint iterationStep, bool entrapBubbles,
        plint numRememberedVolumes);
    void computeNewBubbles(
        std::set<plint> &oldIDs, std::set<plint> &newIDs, std::vector<double> const &bubbleVolume,
        T newBubbleVolumeCorrection, std::map<plint, BubbleInfo3D> &newBubbles,
        std::map<plint, plint> &newToFinal, bool entrapBubbles, plint numRememberedVolumes);
    void updateBubbleLog(
        BubbleTransition3D &bubbleTransition, std::vector<double> const &bubbleVolume,
        plint iterationStep, std::map<plint, BubbleInfo3D> &newBubbles,
        std::map<plint, plint> &newToFinal);
    static void computeBubbleTransitions(
        std::vector<std::vector<plint> > const &newToAllOldMap,
        std::map<plint, std::vector<plint> > const &oldToAllNewMap,
        std::vector<BubbleTransition3D> &bubbleTransitions);
    static T getCutoffVolume()
    {
        return 2.0;
    }

private:
    BubbleHistory3D(BubbleHistory3D<T> const &rhs)
    {
        PLB_ASSERT(false);
    }
    BubbleHistory3D<T> &operator=(BubbleHistory3D<T> const &rhs)
    {
        PLB_ASSERT(false);
        return *this;
    }

private:
    MultiContainerBlock3D *bubbleAnalysisContainer, *bubbleCorrelationContainer,
        *bubbleRemapContainer;
    BubbleMPIdata mpiData;
    MultiScalarField3D<plint> *oldTagMatrix;
    std::map<plint, BubbleInfo3D> bubbles;
    plint nextBubbleID;
    //    <Iteration>   <IDs of created bubbles>,<IDs of vanished bubbles>
    std::map<plint, std::pair<std::vector<plint>, std::vector<plint> > > timeHistory;
    std::vector<FullBubbleRecord> fullBubbleRecord;
    /// IMPORTANT ///
    // Checkpointing needs to be implemented for this code. To checkpoint, you need
    // to save
    // - bubbles
    // - fullBubbleRecord
    // - timeHistory (optional)
    // - oldTagMatrix
    // - nextBubbleID
};

struct BubbleCorrelationData3D : public ContainerBlockData {
    virtual BubbleCorrelationData3D *clone() const
    {
        return new BubbleCorrelationData3D(*this);
    }
    // Maps new bubble tags to old bubble tags at successive iterations.
    // This is needed because a new tag may map to several old tags, due
    // to a bubble merging process.
    std::vector<plint> newToOldMap0, newToOldMap1;
};

class BubbleInfo3D {
public:
    BubbleInfo3D() : referenceVolume(0.), currentVolume(0.), frozen(false)
    {
        numRememberedVolumes = 0;
        rememberedVolumes.clear();
        invNumRememberedVolumes = 0.0;
        numRegisteredVolumes = 0;
        meanVolume = 0.0;
    }
    BubbleInfo3D(double volume, plint numRememberedVolumes_) :
        referenceVolume(volume),
        currentVolume(volume),
        frozen(false),
        numRememberedVolumes(numRememberedVolumes_)
    {
        PLB_ASSERT(numRememberedVolumes > 0);

        if (numRememberedVolumes != 1) {
            rememberedVolumes.resize(numRememberedVolumes, referenceVolume);
            invNumRememberedVolumes = 1.0 / numRememberedVolumes;
            numRegisteredVolumes = 0;
            meanVolume = referenceVolume;
        } else {
            rememberedVolumes.clear();
            invNumRememberedVolumes = 0.0;
            numRegisteredVolumes = 0;
            meanVolume = referenceVolume;
        }
    }
    void freeze()
    {
        frozen = true;
    }
    void setVolume(double newVolume)
    {
        if (!frozen) {
            currentVolume = newVolume;

            if (numRememberedVolumes != 1) {
                plint i = numRegisteredVolumes % numRememberedVolumes;
                numRegisteredVolumes++;
                meanVolume += invNumRememberedVolumes * (currentVolume - rememberedVolumes[i]);
                rememberedVolumes[i] = currentVolume;
            } else {
                meanVolume = currentVolume;
            }
        }
    }
    double getVolumeRatio() const
    {
        static const double epsilon = std::numeric_limits<double>::epsilon() * 1.e4;
        if (std::fabs(meanVolume) > epsilon) {
            return referenceVolume / meanVolume;
        } else {
            return 1.0;
        }
    }
    double getReferenceVolume() const
    {
        return referenceVolume;
    }
    double getVolume() const
    {
        return currentVolume;
    }
    bool isFrozen() const
    {
        return frozen;
    }

private:
    double referenceVolume;
    double currentVolume;
    bool frozen;

    // The following variables relate to the hysteresis mechanism.
    plint numRememberedVolumes;
    std::vector<double> rememberedVolumes;
    double invNumRememberedVolumes;
    plint numRegisteredVolumes;
    double meanVolume;
};

struct BubbleTransition3D {
    bool empty() const
    {
        return oldIDs.empty() && newIDs.empty();
    }
    std::string description() const
    {
        std::stringstream sstream;
        if (oldIDs.empty() && newIDs.empty()) {
        } else if (oldIDs.empty()) {
            PLB_ASSERT(newIDs.size() == 1);
            sstream << "Bubble " << *newIDs.begin() << " created";
        } else if (newIDs.empty()) {
            PLB_ASSERT(oldIDs.size() == 1);
            sstream << "Bubble " << *oldIDs.begin() << " vanished";
        } else if (oldIDs.size() == 1) {
            if (newIDs.size() == 1) {
                sstream << "Straight transition from ID " << *oldIDs.begin() << " to ID "
                        << *newIDs.begin();
            } else {
                sstream << "Splitting from ID " << *oldIDs.begin() << " to IDs";
                std::set<plint>::const_iterator it = newIDs.begin();
                for (; it != newIDs.end(); ++it) {
                    sstream << " " << *it;
                }
            }
        } else if (newIDs.size() == 1) {
            sstream << "Merging IDs";
            std::set<plint>::const_iterator it = oldIDs.begin();
            for (; it != oldIDs.end(); ++it) {
                sstream << " " << *it;
            }
            sstream << " into ID " << *newIDs.begin();
        } else {
            sstream << "Transition from IDs";
            std::set<plint>::const_iterator it1 = oldIDs.begin();
            for (; it1 != oldIDs.end(); ++it1) {
                sstream << " " << *it1;
            }
            sstream << " into IDs ";
            std::set<plint>::const_iterator it2 = newIDs.begin();
            for (; it2 != newIDs.end(); ++it2) {
                sstream << " " << *it2;
            }
        }
        return sstream.str();
    }
    std::set<plint> oldIDs, newIDs;
};

struct FullBubbleRecord {
    FullBubbleRecord(double initialVolume_, plint beginIteration_) :
        initialVolume(initialVolume_),
        finalVolume(0.),
        beginIteration(beginIteration_),
        endIteration(beginIteration_),
        frozen(false)
    { }

    std::string description(plint ID) const
    {
        std::stringstream sstream;
        sstream << "Bubble " << ID << ". At it. " << beginIteration << ", Vol. " << initialVolume
                << ": " << beginTransition.description() << ".";
        if (endIteration > beginIteration) {
            sstream << " Removed at it. " << endIteration << ", Vol. " << finalVolume << ": "
                    << endTransition.description() << ".";
        }
        sstream << " Frozen bubble: " << (frozen ? "Yes." : "No.");
        return sstream.str();
    }

    double initialVolume, finalVolume;
    plint beginIteration, endIteration;
    BubbleTransition3D beginTransition, endTransition;
    bool frozen;
};

template <typename T>
class CorrelateBubbleIds3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual CorrelateBubbleIds3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;  // tags before.
        modified[1] = modif::nothing;  // tags after.
        modified[2] = modif::nothing;  // data.
    }
};

template <typename T>
class UpdateBubblePressure3D : public BoxProcessingFunctional3D_SS<plint, T> {
public:
    UpdateBubblePressure3D(
        std::map<plint, BubbleInfo3D> const &bubbles_, T rho0_, T alpha_, T beta_, T gamma_);
    virtual void process(Box3D domain, ScalarField3D<plint> &tags, ScalarField3D<T> &density);
    virtual UpdateBubblePressure3D *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;          // tags.
        modified[1] = modif::staticVariables;  // density.
    }

private:
    static T getCutoffDensityRatio()
    {
        return 2.0;
    }
    // static T getCutoffDensityRatio() { return 1.01; }
private:
    std::map<plint, BubbleInfo3D> bubbles;
    T rho0;
    T alpha;
    T beta;
    T gamma;
};

}  // namespace plb

#endif  // BUBBLE_HISTORY_3D_H
