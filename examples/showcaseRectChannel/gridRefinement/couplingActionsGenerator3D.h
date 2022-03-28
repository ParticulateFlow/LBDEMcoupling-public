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
 * Coupling between grids of different refinement level -- header file.
 */

#ifndef COUPLING_ACTIONS_GENERATOR_3D_H
#define COUPLING_ACTIONS_GENERATOR_3D_H

#include <vector>

#include "core/globalDefs.h"
#include "gridRefinement/couplingInterfaceGenerator3D.h"
#include "multiBlock/coupling3D.h"
#include "multiBlock/group3D.h"

namespace plb {

// Contains all the necessary information for a complete grid refined lattice.
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
class MultiLevelActions3D : public MultiLevel3D {
public:
    // The dynamics dyn must be right for the coarsest instantiated level.
    // MultiLevelActions3D(OctreeGridStructure& ogs_, Dynamics<T,Descriptor> *dyn_, plint order_,
    //                    plint overlapWidth_ = 1, bool filterAll = true);
    MultiLevelActions3D(
        OctreeGridStructure &ogs_, Dynamics<T, Descriptor> *dyn_, plint order_,
        std::vector<Group3D *> &groups_, plint overlapWidth_ = 1, bool filterAll = true);
    MultiLevelActions3D(MultiLevelActions3D<T, Descriptor, Engine> const &rhs);
    MultiLevelActions3D<T, Descriptor, Engine> &operator=(
        MultiLevelActions3D<T, Descriptor, Engine> const &rhs);
    void swap(MultiLevelActions3D<T, Descriptor, Engine> &rhs);
    ~MultiLevelActions3D();
    void registerBlocks(Actions3D &actions, plint level);

    void initialize();

    void writeInterfaces(double dx, const Array<double, 3> &pos) const
    {
        for (plint iL = 0; iL < (plint)interfaces.size(); ++iL) {
            interfaces[iL]->writeInterfaces(dx, pos);
        }
    }

    virtual MultiBlockLattice3D<T, Descriptor> const &getLevel(plint iL) const
    {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return groups[iL]->template getLattice<T, Descriptor>("lattice");
    }

    virtual MultiBlockLattice3D<T, Descriptor> &getLevel(plint iL)
    {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return groups[iL]->template getLattice<T, Descriptor>("lattice");
    }

    plint getNumLevels() const
    {
        return ogs.getNumLevels();
    }

    plint const &getOrder() const
    {
        return order;
    }

    OctreeGridStructure const &getOgs() const
    {
        return ogs;
    }

    void initializeTensorFields();
    /* std::vector<Actions3D*> actions = generateActions();
     * Actions3D* actions_iL = actions[iL] holds the actions for level iL.
     * Actions3D& step0 = actions_iL->getActions(0);
     * The steps are:
     * STEP0: collide-and-stream, and internal processors.
     * STEP1: Decompose->T0; TimeInterpolate->T12
     * STEP2: RECURSE
     * STEP3: SpaceInterpolate:T12->DecFine; Recompose:DecFine->Lattice.level+1
     * STEP4: RECURSE
     * STEP5: SpaceInterpolate:T0->DecFine; etc...
     *
     * You can add actions to any step. If you add them to a recursion step,
     * your actions will always be executed BEFORE recursion.
     */
    std::vector<Actions3D *> generateActions();
    void addDefaultCouplings(std::vector<Actions3D *> &actions);
    void addCollideAndStream(std::vector<Actions3D *> &actions);
    void addInternalProcessors(std::vector<Actions3D *> &actions);

private:
    void integrateCopyAndSpatialInterpolation(
        MultiNTensorField3D<T> &cTensor, MultiNTensorField3D<T> &fTensor,
        CouplingInterfaces3D &coupling, plint numExec);
    Actions3D *createCopyAndSpatialInterpolation(
        MultiNTensorField3D<T> &cTensor, MultiNTensorField3D<T> &fTensor,
        CouplingInterfaces3D &coupling);

private:
    OctreeGridStructure ogs;
    Dynamics<T, Descriptor> *dyn;
    plint order;
    plint overlapWidth;

private:
    std::vector<Group3D *> groups;
    std::vector<CouplingInterfaces3D *> interfaces;
    bool filterAll;
};

class MultiGridIterator3D {
public:
    MultiGridIterator3D(int numLevels_, std::vector<Actions3D *> &actions_);
    int step();
    void iterate(int whichLevel = 0);

private:
    int numLevels;
    std::vector<Actions3D *> &actions;
    std::vector<int> position;
    int level;
};

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
std::vector<Group3D *> generateLattices(
    Dynamics<T, Descriptor> *dyn, OctreeGridStructure const &ogs, plint envelopeWidth = 1);

template <typename T, template <typename U> class Descriptor>
void generateHelperBlocks(
    std::vector<Group3D *> &groups, OctreeGridStructure const &ogs, plint order);

template <typename T, template <typename U> class Descriptor>
void generateHelperBlocks(
    std::vector<Group3D *> &groups, plint level, OctreeGridStructure const &ogs, plint order);

std::vector<Box3D> optimizeCoarseBoxes(
    std::vector<boxLogic::DirectedPlane> coarseExt, OctreeGridStructure const &ogs, plint level);

std::vector<Box3D> optimizeFineBoxes(
    std::vector<Box3D> fineBoxes, OctreeGridStructure const &ogs, plint fineLevel);

}  // namespace plb

#endif  // COUPLING_ACTIONS_GENERATOR_3D_H
