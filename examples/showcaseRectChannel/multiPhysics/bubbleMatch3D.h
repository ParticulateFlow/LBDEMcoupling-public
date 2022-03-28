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

#ifndef BUBBLE_MATCH_3D_H
#define BUBBLE_MATCH_3D_H

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "core/globalDefs.h"
#include "offLattice/makeSparse3D.h"
#include "parallelism/mpiManager.h"

namespace plb {

class BubbleMPIdata {
public:
    BubbleMPIdata(MultiBlock3D &block);
    std::vector<plint> const &getLocalIds() const;

private:
    void computeLocalIds(MultiBlock3D &block);

private:
    std::vector<plint> localIds;
};

class BubbleMatch3D {
public:
    BubbleMatch3D(MultiBlock3D &templ);
    ~BubbleMatch3D();
    template <typename T>
    void execute(MultiScalarField3D<int> &flag, MultiScalarField3D<T> &volumeFraction);
    MultiScalarField3D<plint> *getTagMatrix()
    {
        return tagMatrix;
    }
    void setTagMatrix(MultiScalarField3D<plint> *newTagMatrix)
    {
        tagMatrix = newTagMatrix;
    }
    std::vector<double> const &getBubbleVolume()
    {
        return bubbleVolume;
    }
    std::vector<Array<double, 3> > const &getBubbleCenter()
    {
        return bubbleCenter;
    }
    pluint numBubbles() const
    {
        return bubbleVolume.size();
    }

private:
    // Re-assign a continuously numbered ID to the detected bubbles.
    pluint countAndTagBubbles();
    // Computes the volumes and centers of all new bubbles.
    template <typename T>
    void bubbleAnalysis(
        MultiScalarField3D<int> &flag, MultiScalarField3D<T> &volumeFraction, pluint numBubbles);
    // Implements all required MPI operations needed to compute the bubble volume and centers,
    // after calling AnalyzeBubbles3D.
    void computeBubbleData(pluint numBubbles);
    // Implements all required MPI operations needed to compute the global IDs of the current
    // bubbbles, after calling CollectBubbleTags3D.
    plint globalBubbleIds();
    // Prepare the bubble container for the next time iteration.
    void resetBubbleContainer();
    // Parallel bucket-fill algorithm to assign a unique ID to every contiguous region.
    void bubbleBucketFill(MultiScalarField3D<int> &flag);

private:
    // TODO: Why is the copy constructor and = doing nothing?
    BubbleMatch3D(BubbleMatch3D const &rhs) : mpiData(rhs.mpiData)
    {
        PLB_ASSERT(false);
    }
    BubbleMatch3D &operator=(BubbleMatch3D const &rhs)
    {
        PLB_ASSERT(false);
        return *this;
    }

private:
    MultiContainerBlock3D *bubbleContainer, *bubbleAnalysisContainer, *bubbleRemapContainer;
    BubbleMPIdata mpiData;
    MultiScalarField3D<plint> *tagMatrix;
    std::vector<double> bubbleVolume;
    std::vector<Array<double, 3> > bubbleCenter;
    static const plint maxNumBubbles = 100000;
};

/**
 * Data for the bubble counter, associated to one block.
 * It holds information that changes during time:
 *  nextCellId: next available ID, if a new cell-type is found.
 *  retagging: a map that relates equivalent IDs, when two domains
 *             are found to be contiguous.
 *  maxNumBubbles: an upper bound for the number of bubbles per block,
 *                 so every block can create a globally unique bubble ID.
 **/
class BubbleCounterData3D : public ContainerBlockData {
public:
    BubbleCounterData3D(plint maxNumBubbles_);
    virtual BubbleCounterData3D *clone() const;
    // The assumption here is that cell 0 is a cell that needs to
    // be tagged ("a bubble cell"). Depending on the tag of neighboring
    // cells, either one of the neighbor's tag or a new tag is assigned.
    // There is a conflict if cell 0 has been previously converted in
    // a way that is incompatible with the neighbor's tags.
    bool convertCell(
        plint &tag0, plint tag1, plint tag2, plint tag3, plint tag4, plint tag5, plint tag6,
        plint tag7, plint tag8, plint tag9, plint tag10, plint tag11, plint tag12, plint tag13,
        plint tag1_, plint tag2_, plint tag3_, plint tag4_, plint tag5_, plint tag6_, plint tag7_,
        plint tag8_, plint tag9_, plint tag10_, plint tag11_, plint tag12_, plint tag13_);
    // The assumption here is that cell 0 is a cell that needs to
    // be tagged ("a bubble cell"). Depending on the tag of neighboring
    // cell 1, either cell 1's tag or a new tag is assigned.
    // There is a conflict if cell 0 has been previously converted in
    // a way that is incompatible with cell 1's tag.
    bool processNeighbor(plint &tag0, plint tag1);
    plint getNextTag();
    void reset();
    // Get the "real" value (after remapping) of a given tag.
    plint convertTag(plint tag) const;
    // Important: it is assumed that oldTag and newTag are not
    // themselves mapped. This means that for arbitrary tags tag1 and
    // tag2 you should not call registerConflict(tag1,tag2), but
    // registerConflict(convertTag(tag1), convertTag(tag2).

    void registerConflict(plint oldTag, plint newTag);

private:
    plint nextCellId;
    std::map<plint, plint> retagging;

private:
    plint maxNumBubbles;
};

class BubbleRemapData3D : public ContainerBlockData {
public:
    BubbleRemapData3D(plint maxNumBubbles_ = 0) : maxNumBubbles(maxNumBubbles_) { }
    virtual BubbleRemapData3D *clone() const;
    std::vector<plint> &getUniqueTags()
    {
        return uniqueTags;
    }
    std::vector<plint> const &getUniqueTags() const
    {
        return uniqueTags;
    }
    std::map<plint, plint> &getTagRemap()
    {
        return tagRemap;
    }
    bool isMyTag(plint tag);

private:
    plint maxNumBubbles;
    std::vector<plint> uniqueTags;
    std::map<plint, plint> tagRemap;
};

struct BubbleAnalysisData3D : public ContainerBlockData {
    virtual BubbleAnalysisData3D *clone() const
    {
        return new BubbleAnalysisData3D(*this);
    }
    std::vector<double> bubbleVolume;
    std::vector<Array<double, 3> > bubbleCenter;
};

class CountBubbleIteration3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    CountBubbleIteration3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual CountBubbleIteration3D *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;  // tags.
        modified[1] = modif::nothing;          // flags.
        modified[2] = modif::nothing;          // data.
    }
    plint getNumConflicts() const;

private:
    plint numConflictsId;
};

template <typename T>
class AnalyzeBubbles3D : public BoxProcessingFunctional3D {
public:
    AnalyzeBubbles3D(pluint numBubbles_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual AnalyzeBubbles3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;  // tags.
        modified[1] = modif::nothing;          // flags.
        modified[2] = modif::nothing;          // data.
        modified[3] = modif::nothing;          // volume fraction.
    }

private:
    pluint numBubbles;
};

// Converts the information about the overall available bubble tags available
// to all processors.
class CollectBubbleTags3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual CollectBubbleTags3D *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;  // tags.
        modified[1] = modif::nothing;  // data.
    }
};

// Assign a new tag to all bubble cells (they must have been uniquely tagged previously).
// The only field in the BubbleCounterData3D which is used here is tagRemap.
class ApplyTagRemap3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual ApplyTagRemap3D *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;  // tags.
        modified[1] = modif::nothing;          // data.
    }
};

}  // namespace plb

#endif  // BUBBLE_MATCH_3D_H
