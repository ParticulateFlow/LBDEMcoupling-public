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

#ifndef COUPLING_ACTIONS_GENERATOR_3D_HH
#define COUPLING_ACTIONS_GENERATOR_3D_HH

#include <map>
#include <vector>

#include "core/globalDefs.h"
#include "gridRefinement/boxLogic3D.hh"
#include "gridRefinement/couplingActionsGenerator3D.h"
#include "gridRefinement/gridRefinementFunctional3D.h"
#include "gridRefinement/gridRefinementFunctional3D.hh"
#include "gridRefinement/octreeGridStructure.h"
#include "gridRefinement/rescaleEngine.h"

namespace plb {

/*
template<typename T, template <typename U> class Descriptor,
template<typename T2, template<typename U2> class Descriptor2> class Engine>
MultiLevelActions3D<T,Descriptor,Engine>::MultiLevelActions3D (
OctreeGridStructure& ogs_, Dynamics<T,Descriptor> *dyn_, plint order_, plint overlapWidth_, bool
filterAll_ ) :   ogs(ogs_), dyn(dyn_), order(order_), overlapWidth(overlapWidth_),
    filterAll(filterAll_)
{
pcout << "Generate levels and coupling-interfaces" << std::endl;
generateLevels();
// Creation of the coupling interfaces vector
for (plint level=0; level < getNumLevels()-1; ++level) {
    interfaces.push_back(new CouplingInterfaces3D(ogs, level+1,overlapWidth));
}
pcout << "initializeTensorFields" << std::endl;
initializeTensorFields();
}
*/

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
MultiLevelActions3D<T, Descriptor, Engine>::MultiLevelActions3D(
    OctreeGridStructure &ogs_, Dynamics<T, Descriptor> *dyn_, plint order_,
    std::vector<Group3D *> &groups_, plint overlapWidth_, bool filterAll_) :
    ogs(ogs_),
    dyn(dyn_),
    order(order_),
    overlapWidth(overlapWidth_),
    groups(groups_),
    filterAll(filterAll_)
{
    // pcout << "Generate coupling-interfaces" << std::endl;
    //  Creation of the coupling interfaces vector
    for (plint level = 0; level < getNumLevels() - 1; ++level) {
        interfaces.push_back(new CouplingInterfaces3D(ogs, level + 1, overlapWidth));
    }
    // pcout << "Initialize tensor-fields" << std::endl;
    initializeTensorFields();
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
MultiLevelActions3D<T, Descriptor, Engine>::MultiLevelActions3D(
    MultiLevelActions3D<T, Descriptor, Engine> const &rhs) :
    ogs(rhs.ogs),
    dyn(rhs.dyn),
    order(rhs.order),
    overlapWidth(rhs.overlapWidth),
    interfaces(rhs.interfaces),
    filterAll(rhs.filterAll)
{ }

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
MultiLevelActions3D<T, Descriptor, Engine> &MultiLevelActions3D<T, Descriptor, Engine>::operator=(
    MultiLevelActions3D<T, Descriptor, Engine> const &rhs)
{
    MultiLevelActions3D<T, Descriptor, Engine>(rhs).swap(*this);
    return *this;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelActions3D<T, Descriptor, Engine>::swap(
    MultiLevelActions3D<T, Descriptor, Engine> &rhs)
{
    std::swap(ogs, rhs.ogs);
    std::swap(dyn, rhs.dyn);
    std::swap(order, rhs.order);
    std::swap(overlapWidth, rhs.overlapWidth);
    std::swap(interfaces, rhs.interfaces);
    std::swap(filterAll, rhs.filterAll);
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
MultiLevelActions3D<T, Descriptor, Engine>::~MultiLevelActions3D()
{
    delete dyn;
    for (pluint iA = 0; iA < interfaces.size(); ++iA) {
        delete interfaces[iA];
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelActions3D<T, Descriptor, Engine>::registerBlocks(Actions3D &actions, plint level)
{
    actions.addBlock(groups[level]->template getLattice<T, Descriptor>("lattice"));      // 0
    actions.addBlock(groups[level]->template getNTensor<T>("decomposed_t0"));            // 1
    actions.addBlock(groups[level]->template getNTensor<T>("decomposed_t1"));            // 2
    actions.addBlock(groups[level]->template getNTensor<T>("decomposed_t12"));           // 3
    actions.addBlock(groups[level + 1]->template getNTensor<T>("decomposed_fine"));      // 4
    actions.addBlock(groups[level + 1]->template getLattice<T, Descriptor>("lattice"));  // 5
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelActions3D<T, Descriptor, Engine>::initialize()
{
    // initialization of the lattice in each grid level
    PLB_ASSERT(getNumLevels() == (plint)groups.size());
    for (pluint iL = 0; iL < getNumLevels(); ++iL) {
        groups[iL]->template getLattice<T, Descriptor>("lattice").initialize();
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelActions3D<T, Descriptor, Engine>::initializeTensorFields()
{
    // initialization of the t0 NTenforfields
    for (plint iL = 0; iL < getNumLevels() - 1; ++iL) {
        std::vector<boxLogic::DirectedPlane> cToFext_cu = interfaces[iL]->getCtoFext_cu();
        for (pluint iA = 0; iA < cToFext_cu.size(); ++iA) {
            applyProcessingFunctional(
                new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>((T)0.5, T(), order),
                cToFext_cu[iA].bb, groups[iL]->template getLattice<T, Descriptor>("lattice"),
                groups[iL]->template getNTensor<T>("decomposed_t0"));
            applyProcessingFunctional(
                new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>((T)0.5, T(), order),
                cToFext_cu[iA].bb, groups[iL]->template getLattice<T, Descriptor>("lattice"),
                groups[iL]->template getNTensor<T>("decomposed_t1"));
        }
        std::vector<Box3D> fToCint_fu = interfaces[iL]->getFtoCint_fu();
        for (pluint iA = 0; iA < fToCint_fu.size(); ++iA) {
            applyProcessingFunctional(
                new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>((T)2, T(), order),
                fToCint_fu[iA], groups[iL + 1]->template getLattice<T, Descriptor>("lattice"),
                groups[iL + 1]->template getNTensor<T>("decomposed_fine"));
        }
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
std::vector<Actions3D *> MultiLevelActions3D<T, Descriptor, Engine>::generateActions()
{
    // initialization of the t0 NTenforfields
    std::vector<Actions3D *> generatedActions(getNumLevels());
    for (plint iL = 0; iL < getNumLevels(); ++iL) {
        Actions3D *level = new Actions3D();
        Actions3D *step0 = new Actions3D();
        if (iL == getNumLevels() - 1) {
            step0->addBlock(groups[iL]->template getLattice<T, Descriptor>("lattice"));
            level->addActions(step0);
        } else {
            Actions3D *step1 = new Actions3D();
            Actions3D *step3 = new Actions3D();
            Actions3D *step5 = new Actions3D();
            registerBlocks(*step0, iL);
            registerBlocks(*step1, iL);
            registerBlocks(*step3, iL);
            registerBlocks(*step5, iL);
            level->addActions(step0);
            level->addActions(step1);
            level->addNOP();  // Step 2
            level->addActions(step3);
            level->addNOP();  // Step 4
            level->addActions(step5);
        }
        generatedActions[iL] = level;
    }
    return generatedActions;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelActions3D<T, Descriptor, Engine>::addCollideAndStream(
    std::vector<Actions3D *> &actions)
{
    PLB_ASSERT((plint)actions.size() == getNumLevels());
    for (plint iL = 0; iL < getNumLevels(); ++iL) {
        Actions3D &step0 = actions[iL]->getActions(0);
        step0.addCollideAndStream<T, Descriptor>(0);
        step0.addAllInternalProcessorsAndComm(0);
        step0.addIncrementTime<T, Descriptor>(0);
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelActions3D<T, Descriptor, Engine>::addInternalProcessors(
    std::vector<Actions3D *> &actions)
{
    PLB_ASSERT((plint)actions.size() == getNumLevels());
    for (pluint iL = 0; iL < getNumLevels(); ++iL) {
        Actions3D &step0 = actions[iL]->getActions(0);
        step0.addAllInternalProcessorsAndComm(0);
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelActions3D<T, Descriptor, Engine>::addDefaultCouplings(
    std::vector<Actions3D *> &actions)
{
    PLB_ASSERT((plint)actions.size() == getNumLevels());
    for (plint iL = 0; iL < getNumLevels() - 1; ++iL) {  // Don't add interfaces to finest level
        /***** STEP 1: Decompose->T0; TimeInterpolate->T12 ***********/
        Actions3D &step1 = actions[iL]->getActions(1);
        step1.addProcessor(
            new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>((T)0.5, T(), order), 0, 1,
            interfaces[iL]->boxCtoFext_cu());
        step1.addCommunication(1, modif::staticVariables);  // decomposedT0
        // Applies to bulk and envelope: no communication needed.
        step1.addProcessor(
            new TemporalInterpolationFunctional3D<T>(), 1, 2, 3, interfaces[iL]->boxCtoFext_cu());
        /***** STEP 2: RECURSE ***********/
        /***** STEP 3: SpaceInterpolate:T12->DecFine; Recompose:DecFine->Lattice.level+1
         * ***********/
        Actions3D &step3 = actions[iL]->getActions(3);
        step3.addActions(createCopyAndSpatialInterpolation(
            groups[iL]->template getNTensor<T>("decomposed_t12"),
            groups[iL + 1]->template getNTensor<T>("decomposed_fine"), *interfaces[iL]));
        step3.addCommunication(4, modif::staticVariables);  // decomposedFine
        // Applies to bulk and envelope: no communication needed.
        step3.addProcessor(
            new RecomposeFunctional3D<T, Descriptor>(order), 5, 4, interfaces[iL]->boxCtoFint_fu());
        /***** STEP 4: RECURSE ***********/
        /***** STEP 5: SpaceInterpolate:T0->DecFine; Recompose:DecFine->Lattice.level+1;
         * Decompose:Lattice.level+1->DecFine; Filter/Copy:DecFine->T12; Recompose:T12->Lattice;
         * Swap: T0<->T1  ***********/
        Actions3D &step5 = actions[iL]->getActions(5);
        step5.addActions(createCopyAndSpatialInterpolation(
            groups[iL]->template getNTensor<T>("decomposed_t0"),
            groups[iL + 1]->template getNTensor<T>("decomposed_fine"), *interfaces[iL]));
        step5.addCommunication(4, modif::staticVariables);  // decomposedFine
        step5.addProcessor(
            new RecomposeFunctional3D<T, Descriptor>(order), 5, 4, interfaces[iL]->boxCtoFint_fu());
        step5.addProcessor(
            new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>((T)2, T(), order), 5, 4,
            interfaces[iL]->getFtoCint_fu());
        step5.addCommunication(4, modif::staticVariables);  // decomposedFine
        // Don't communicate after this step, because we communicate for T12(block 3) after
        // CopyFunctional.
        step5.addProcessor(
            new CopyAndFilterNonEquilibriumFunctional3D<T, Descriptor>(filterAll), 3, 4,
            interfaces[iL]->boxFtoCint_cu());
        step5.addProcessor(
            new CopyFunctional3D<T, Descriptor>(), 3, 4, interfaces[iL]->boxFtoCboundaryEdges_cu());
        step5.addCommunication(3, modif::staticVariables);
        // No communication needed: applies to bulk and envelope.
        step5.addProcessor(
            new RecomposeFunctional3D<T, Descriptor>(order), 0, 3, interfaces[iL]->boxFtoCint_cu());
        // No communication needed: applies to bulk and envelope.
        step5.addProcessor(
            new RecomposeFunctional3D<T, Descriptor>(order), 0, 3,
            interfaces[iL]->boxFtoCboundaryEdges_cu());
        // No communication needed: applies to bulk and envelope.
        step5.addProcessor(
            new SwapValuesBulkAndEnvelope3D_N<T>(), 1, 2,
            groups[iL]->template getNTensor<T>("decomposed_t0").getBoundingBox());
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelActions3D<T, Descriptor, Engine>::integrateCopyAndSpatialInterpolation(
    MultiNTensorField3D<T> &cTensor, MultiNTensorField3D<T> &fTensor,
    CouplingInterfaces3D &coupling, plint numExec)
{
    // pcout << "Bulk Planes = " << std::endl;
    for (pluint iA = 0; iA < coupling.getCtoFplanes_cu().size(); ++iA) {
        Box3D bb = coupling.getCtoFplanes_cu()[iA].bb;
        // pcout << "iA = " << iA << " " << coupling.getCtoFplanes_cu()[iA].toStr() << std::endl;
        plint dir = coupling.getCtoFplanes_cu()[iA].direction;
        integrateProcessingFunctional(
            new CopyAndSpatialInterpolationPlaneFunctional3D<T>(dir), bb, cTensor, fTensor,
            numExec);
    }

    // pcout << "Edges for bulk = " << std::endl;
    for (pluint iA = 0; iA < coupling.getCtoFbulkEdges_cu().size(); ++iA) {
        Box3D bb = coupling.getCtoFbulkEdges_cu()[iA].bb;
        // pcout << "iA = " << iA << " " << printBox(bb) << std::endl;
        plint planeDir = coupling.getCtoFbulkEdges_cu()[iA].planeDir;
        plint dir1 = coupling.getCtoFbulkEdges_cu()[iA].dir1;
        plint dir2 = coupling.getCtoFbulkEdges_cu()[iA].dir2;
        integrateProcessingFunctional(
            new CopyAndSpatialInterpolationEdgeFunctional3D<T>(planeDir, dir1, dir2), bb, cTensor,
            fTensor, numExec);
    }

    // pcout << "Corners for bulk = " << std::endl;
    for (pluint iA = 0; iA < coupling.getCtoFbulkCorners_cu().size(); ++iA) {
        Box3D bb = coupling.getCtoFbulkCorners_cu()[iA].bb;
        // pcout << "iA = " << iA << " " << printBox(bb) << std::endl;
        plint planeDir = coupling.getCtoFbulkCorners_cu()[iA].planeDir;
        plint dir1 = coupling.getCtoFbulkCorners_cu()[iA].dir1;
        plint dir2 = coupling.getCtoFbulkCorners_cu()[iA].dir2;
        integrateProcessingFunctional(
            new CopyAndSpatialInterpolationCornerFunctional3D<T>(planeDir, dir1, dir2), bb, cTensor,
            fTensor, numExec);
    }

    // pcout << "Edges for BCs = " << std::endl;
    for (pluint iA = 0; iA < coupling.getCtoFboundaryEdges_cu().size(); ++iA) {
        Box3D bb = coupling.getCtoFboundaryEdges_cu()[iA].bb;
        // pcout << "iA = " << iA << " " << printBox(bb) << std::endl;
        plint planeDir = coupling.getCtoFboundaryEdges_cu()[iA].planeDir;
        plint dir1 = coupling.getCtoFboundaryEdges_cu()[iA].dir1;
        plint dir2 = coupling.getCtoFboundaryEdges_cu()[iA].dir2;
        integrateProcessingFunctional(
            new CopyAndSpatialInterpolationBoundaryEdgeFunctional3D<T>(planeDir, dir1, dir2), bb,
            cTensor, fTensor, numExec);
    }

    // pcout << "Corners for BCs = " << std::endl;
    for (pluint iA = 0; iA < coupling.getCtoFboundaryCorners_cu().size(); ++iA) {
        Box3D bb = coupling.getCtoFboundaryCorners_cu()[iA].bb;
        // pcout << "iA = " << iA << " " << printBox(bb) << std::endl;
        plint planeDir = coupling.getCtoFboundaryCorners_cu()[iA].planeDir;
        plint dir1 = coupling.getCtoFboundaryCorners_cu()[iA].dir1;
        plint dir2 = coupling.getCtoFboundaryCorners_cu()[iA].dir2;

        integrateProcessingFunctional(
            new CopyAndSpatialInterpolationBoundaryCornerFunctional3D<T>(planeDir, dir1, dir2), bb,
            cTensor, fTensor, numExec);
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
Actions3D *MultiLevelActions3D<T, Descriptor, Engine>::createCopyAndSpatialInterpolation(
    MultiNTensorField3D<T> &cTensor, MultiNTensorField3D<T> &fTensor,
    CouplingInterfaces3D &coupling)
{
    Actions3D *newActions = new Actions3D();
    newActions->addBlock(cTensor);
    newActions->addBlock(fTensor);
    for (pluint iA = 0; iA < coupling.getCtoFplanes_cu().size(); ++iA) {
        Box3D bb = coupling.getCtoFplanes_cu()[iA].bb;
        plint dir = coupling.getCtoFplanes_cu()[iA].direction;
        newActions->addProcessor(
            new CopyAndSpatialInterpolationPlaneFunctional3D<T>(dir), 0, 1, bb);
    }

    for (pluint iA = 0; iA < coupling.getCtoFbulkEdges_cu().size(); ++iA) {
        Box3D bb = coupling.getCtoFbulkEdges_cu()[iA].bb;
        plint planeDir = coupling.getCtoFbulkEdges_cu()[iA].planeDir;
        plint dir1 = coupling.getCtoFbulkEdges_cu()[iA].dir1;
        plint dir2 = coupling.getCtoFbulkEdges_cu()[iA].dir2;
        newActions->addProcessor(
            new CopyAndSpatialInterpolationEdgeFunctional3D<T>(planeDir, dir1, dir2), 0, 1, bb);
    }

    for (pluint iA = 0; iA < coupling.getCtoFbulkCorners_cu().size(); ++iA) {
        Box3D bb = coupling.getCtoFbulkCorners_cu()[iA].bb;
        plint planeDir = coupling.getCtoFbulkCorners_cu()[iA].planeDir;
        plint dir1 = coupling.getCtoFbulkCorners_cu()[iA].dir1;
        plint dir2 = coupling.getCtoFbulkCorners_cu()[iA].dir2;
        newActions->addProcessor(
            new CopyAndSpatialInterpolationCornerFunctional3D<T>(planeDir, dir1, dir2), 0, 1, bb);
    }

    for (pluint iA = 0; iA < coupling.getCtoFboundaryEdges_cu().size(); ++iA) {
        Box3D bb = coupling.getCtoFboundaryEdges_cu()[iA].bb;
        plint planeDir = coupling.getCtoFboundaryEdges_cu()[iA].planeDir;
        plint dir1 = coupling.getCtoFboundaryEdges_cu()[iA].dir1;
        plint dir2 = coupling.getCtoFboundaryEdges_cu()[iA].dir2;
        newActions->addProcessor(
            new CopyAndSpatialInterpolationBoundaryEdgeFunctional3D<T>(planeDir, dir1, dir2), 0, 1,
            bb);
    }

    for (pluint iA = 0; iA < coupling.getCtoFboundaryCorners_cu().size(); ++iA) {
        Box3D bb = coupling.getCtoFboundaryCorners_cu()[iA].bb;
        plint planeDir = coupling.getCtoFboundaryCorners_cu()[iA].planeDir;
        plint dir1 = coupling.getCtoFboundaryCorners_cu()[iA].dir1;
        plint dir2 = coupling.getCtoFboundaryCorners_cu()[iA].dir2;
        newActions->addProcessor(
            new CopyAndSpatialInterpolationBoundaryCornerFunctional3D<T>(planeDir, dir1, dir2), 0,
            1, bb);
    }
    return newActions;
}
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
std::vector<Group3D *> generateLattices(
    Dynamics<T, Descriptor> *dyn, OctreeGridStructure const &ogs, plint envelopeWidth)
{
    // Computation of the relaxation time vector used for the dynamics.
    std::vector<Array<T, Descriptor<T>::q> > relaxationFrequencies(ogs.getNumLevels());
    Array<T, Descriptor<T>::q> relFreq = dyn->getRelaxationFrequencies();
    Engine<T, Descriptor> engine;

    for (plint level = 0; level < ogs.getNumLevels(); ++level) {
        relaxationFrequencies[level] = relFreq;
        // going from coarse to fine => xDt = 0.5, xInvDx=0 => is useless and
        Array<T, Descriptor<T>::q> relFreqTmp = engine.computeRescaledRelFreq(relFreq, (T)0.5);
        relFreq = relFreqTmp;
    }

    std::vector<Group3D *> lattices;
    for (plint iL = 0; iL < ogs.getNumLevels(); ++iL) {
        dyn->setRelaxationFrequencies(relaxationFrequencies[iL]);
        MultiBlockLattice3D<T, Descriptor> *lattice =
            generateMultiBlockLattice<T, Descriptor>(
                ogs.getMultiBlockManagement(iL, envelopeWidth), dyn->clone())
                .release();
        lattice->periodicity().toggleAll(false);
        defineDynamics(*lattice, lattice->getBoundingBox(), dyn->clone());
        lattices.push_back(new Group3D(lattice, "lattice"));
    }
    delete dyn;
    return lattices;
}

template <typename T, template <typename U> class Descriptor>
void generateHelperBlocks(
    std::vector<Group3D *> &groups, OctreeGridStructure const &ogs, plint order)
{
    for (plint level = 0; level < ogs.getNumLevels() - 1; ++level) {
        generateHelperBlocks<T, Descriptor>(groups, level, ogs, order);
    }
}

template <typename T, template <typename U> class Descriptor>
void generateHelperBlocks(
    std::vector<Group3D *> &groups, plint level, OctreeGridStructure const &ogs, plint order)
{
    PLB_ASSERT(level < ogs.getNumLevels() - 1);

    MultiBlockLattice3D<T, Descriptor> &lattice =
        groups[level]->getLattice<T, Descriptor>("lattice");
    plint numVariables = lattice.getBackgroundDynamics().numDecomposedVariables(order);

    CouplingInterfaces3D interfaces(ogs, level + 1, 1);
    std::vector<Box3D> coarseBoxes = optimizeCoarseBoxes(interfaces.getCtoFext_cu(), ogs, level);
    std::vector<Box3D> fineBoxes = optimizeFineBoxes(interfaces.getFtoCint_fu(), ogs, level + 1);

    plint envCoarse = 2;
    plint envFine = 1;
    bool crop = true;
    groups[level]->generateNTensor<T>(
        "decomposed_t0", numVariables, coarseBoxes, envCoarse, level, crop);
    groups[level]->generateNTensor<T>(
        "decomposed_t12", numVariables, coarseBoxes, envCoarse, level, crop);
    groups[level]->generateNTensor<T>(
        "decomposed_t1", numVariables, coarseBoxes, envCoarse, level, crop);
    groups[level + 1]->generateNTensor<T>(
        "decomposed_fine", numVariables, fineBoxes, envFine, level + 1, crop);
}

}  // namespace plb

#endif  // COUPLING_INTERFACE_GENERATOR_3D_HH
