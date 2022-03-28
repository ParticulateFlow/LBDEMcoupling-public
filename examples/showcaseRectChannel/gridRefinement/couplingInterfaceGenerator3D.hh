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

#ifndef COUPLING_INTERFACE_GENERATOR_3D_HH
#define COUPLING_INTERFACE_GENERATOR_3D_HH

#include <map>
#include <vector>

#include "core/globalDefs.h"
#include "gridRefinement/boxLogic3D.hh"
#include "gridRefinement/couplingInterfaceGenerator3D.h"
#include "gridRefinement/gridRefinementFunctional3D.h"
#include "gridRefinement/gridRefinementFunctional3D.hh"
#include "gridRefinement/octreeGridStructure.h"
#include "gridRefinement/rescaleEngine.h"

namespace plb {

// ======================================================================= //
// ====================GridLevelContainer3D=============================== //
// ======================================================================= //

template <typename T, template <typename U> class Descriptor>
GridLevelContainer3D<T, Descriptor>::GridLevelContainer3D(
    std::vector<Group3D *> &groups, plint order_, plint level_) :
    dyn(0), order(order_), level(level_), parity(0), owned(false)
{
    numVariables = -1;
    PLB_ASSERT(level < (plint)groups.size());
    lattice = &groups[level]->template getLattice<T, Descriptor>("lattice");
    if (level == (plint)groups.size() - 1) {
        decomposed_t0 = decomposed_t12 = decomposed_t1 = decomposed_fine = 0;
    } else {
        decomposed_t0 = &groups[level]->template getNTensor<T>("decomposed_t0");
        decomposed_t12 = &groups[level]->template getNTensor<T>("decomposed_t12");
        decomposed_t1 = &groups[level]->template getNTensor<T>("decomposed_t1");
        decomposed_fine = &groups[level + 1]->template getNTensor<T>("decomposed_fine");
    }
}

template <typename T, template <typename U> class Descriptor>
GridLevelContainer3D<T, Descriptor>::GridLevelContainer3D(
    Dynamics<T, Descriptor> *dyn_, const OctreeGridStructure &ogs, plint order_, plint level_,
    bool xPeriodic_, bool yPeriodic_, bool zPeriodic_) :
    dyn(dyn_), order(order_), level(level_), parity(0), owned(true)
{
    numVariables = dyn->numDecomposedVariables(order);

    lattice = new MultiBlockLattice3D<T, Descriptor>(
        ogs.getMultiBlockManagement(level, 1), defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(), dyn->clone());
    defineDynamics(
        *lattice, lattice->getBoundingBox(),
        dyn->clone());  // TODO: Check if really needed, or should be done a posteriori.
    lattice->periodicity().toggle(0, xPeriodic_);
    lattice->periodicity().toggle(1, yPeriodic_);
    lattice->periodicity().toggle(2, zPeriodic_);
    dataProcessors = new MultiContainerBlock3D(
        lattice->getMultiBlockManagement(), defaultMultiBlockPolicy3D().getCombinedStatistics());

    if (level < ogs.getNumLevels() - 1) {
        plint nextLevel = level + 1;

        CouplingInterfaces3D coupling(ogs, nextLevel, 1);
        std::vector<boxLogic::DirectedPlane> coarseExt = coupling.getCtoFext_cu();
        std::vector<Box3D> coarseBoxes;
        for (std::vector<boxLogic::DirectedPlane>::iterator it = coarseExt.begin();
             it != coarseExt.end(); ++it)
        {
            // Array<plint,2> ary((it->orientation == -1) ? 1 : 0, (it->orientation == 1) ? 1 : 0);
            Array<plint, 6> ary;
            ary[0] = -1;
            ary[1] = 1;
            ary[2] = -1;
            ary[3] = 1;
            ary[4] = -1;
            ary[5] = 1;

            if (it->orientation == 1) {
                ary[2 * it->direction + 1] = 0;
            } else {
                ary[2 * it->direction] = 0;
            }

            // ary[2 * it->direction   ] = (it->orientation == -1) ?  1 : 0;
            // ary[2 * it->direction+1 ] = (it->orientation ==  1) ? -1 : 0;
            // ary[2 * it->direction   ] = -1;
            // ary[2 * it->direction+1 ] = +1;

            coarseBoxes.push_back(it->bb.enlarge(ary));
        }

        // coarseBoxes = boxLogic::removeOverlapsOfBoxes(coarseBoxes);
        // coarseBoxes = boxLogic::merge(coarseBoxes);

        std::vector<Box3D> fineBoxes = coupling.getFtoCint_fu();
        for (std::vector<Box3D>::iterator it = fineBoxes.begin(); it != fineBoxes.end(); ++it) {
            *it = it->enlarge(1);
        }

        // fineBoxes = boxLogic::removeOverlapsOfBoxes(fineBoxes);
        // fineBoxes = boxLogic::merge(fineBoxes);

        Box3D bulk;
        plint levelTmp = 0;
        plint processId = 0;
        std::vector<Box3D> intDom_cu, intDom_fu;
        std::vector<plint> blockIdsPerLevel = ogs.getBlockIdsAtLevel(level, false);
        for (plint iA = 0; iA < (plint)blockIdsPerLevel.size(); ++iA) {
            if (boxLogic::hasFinerNeighbor(ogs, blockIdsPerLevel[iA])) {
                ogs.getBlock(blockIdsPerLevel[iA], bulk, levelTmp, processId);
                intDom_cu.push_back(bulk);
            }
        }
        blockIdsPerLevel = ogs.getOverlapBlockIdsAtLevel(level);
        for (plint iA = 0; iA < (plint)blockIdsPerLevel.size(); ++iA) {
            ogs.getBlock(blockIdsPerLevel[iA], bulk, levelTmp, processId);
            intDom_cu.push_back(bulk);
        }

        // intDom_cu = boxLogic::merge(intDom_cu);

        std::vector<Box3D> finalCoarseBoxes;
        for (plint iA = 0; iA < (plint)intDom_cu.size(); ++iA) {
            std::vector<Box3D> tmpBoxes;
            for (plint iB = 0; iB < (plint)coarseBoxes.size(); ++iB) {
                Box3D newBox;
                bool doesIntersect = intersect(intDom_cu[iA], coarseBoxes[iB], newBox);
                if (doesIntersect) {
                    tmpBoxes.push_back(newBox);
                }
            }
            if (tmpBoxes.size() > 0) {
                tmpBoxes = boxLogic::removeOverlapsOfBoxes(tmpBoxes);
                tmpBoxes = boxLogic::merge(tmpBoxes);
            }
            if (tmpBoxes.size() > 1) {
                Box3D bb = findBoundingBox(tmpBoxes);
                tmpBoxes.clear();
                tmpBoxes.push_back(bb);
            }
            finalCoarseBoxes.insert(finalCoarseBoxes.begin(), tmpBoxes.begin(), tmpBoxes.end());
        }
        finalCoarseBoxes = boxLogic::merge(finalCoarseBoxes);

        blockIdsPerLevel = ogs.getBlockIdsAtLevel(nextLevel, false);
        for (plint iA = 0; iA < (plint)blockIdsPerLevel.size(); ++iA) {
            if (boxLogic::hasCoarserNeighbor(ogs, blockIdsPerLevel[iA])) {
                ogs.getBlock(blockIdsPerLevel[iA], bulk, levelTmp, processId);
                intDom_fu.push_back(bulk);
            }
        }

        // intDom_fu = boxLogic::merge(intDom_fu);

        std::vector<Box3D> finalFineBoxes;
        for (plint iA = 0; iA < (plint)intDom_fu.size(); ++iA) {
            std::vector<Box3D> tmpBoxes;
            for (plint iB = 0; iB < (plint)fineBoxes.size(); ++iB) {
                Box3D newBox;
                bool doesIntersect = intersect(intDom_fu[iA], fineBoxes[iB], newBox);
                if (doesIntersect) {
                    // pcout << "intersection" << std::endl;
                    tmpBoxes.push_back(newBox);
                }
            }
            if (tmpBoxes.size() > 0) {
                tmpBoxes = boxLogic::removeOverlapsOfBoxes(tmpBoxes);
                tmpBoxes = boxLogic::merge(tmpBoxes);
            }
            if (tmpBoxes.size() > 1) {
                Box3D bb = findBoundingBox(tmpBoxes);
                tmpBoxes.clear();
                tmpBoxes.push_back(bb);
            }
            finalFineBoxes.insert(finalFineBoxes.begin(), tmpBoxes.begin(), tmpBoxes.end());
        }
        finalFineBoxes = boxLogic::merge(finalFineBoxes);

        // TODO: make multiblock management consistent with align.
        MultiBlockManagement3D managementCoarse = ogs.getMultiBlockManagement(level, 2);

        // plint envelopeFine = (level == ogs.getNumLevels()-2) ? 1 : 2;
        plint envelopeFine = 1;
        MultiBlockManagement3D managementFine =
            ogs.getMultiBlockManagement(level + 1, envelopeFine);

        MultiBlockManagement3D coarseManagement =
            align(finalCoarseBoxes, managementCoarse, 2, level, true);
        MultiBlockManagement3D fineManagement =
            align(finalFineBoxes, managementFine, envelopeFine, level + 1, true);

        decomposed_t0 = new MultiNTensorField3D<T>(
            numVariables, coarseManagement, defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

        decomposed_t12 = new MultiNTensorField3D<T>(
            numVariables, coarseManagement, defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

        decomposed_t1 = new MultiNTensorField3D<T>(
            numVariables, coarseManagement, defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

        decomposed_fine = new MultiNTensorField3D<T>(
            numVariables, fineManagement, defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

        // decomposed_t0  = new MultiNTensorField3D<T>(numVariables,
        //                                          managementCoarse,
        //                                          defaultMultiBlockPolicy3D().getBlockCommunicator(),
        //                                          defaultMultiBlockPolicy3D().getCombinedStatistics(),
        //                                          defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

        // decomposed_t12 = new MultiNTensorField3D<T>(numVariables,
        //                                       managementCoarse,
        //                                       defaultMultiBlockPolicy3D().getBlockCommunicator(),
        //                                       defaultMultiBlockPolicy3D().getCombinedStatistics(),
        //                                       defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

        // decomposed_t1  = new MultiNTensorField3D<T>(numVariables,
        //                                       managementCoarse,
        //                                       defaultMultiBlockPolicy3D().getBlockCommunicator(),
        //                                       defaultMultiBlockPolicy3D().getCombinedStatistics(),
        //                                       defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

        // decomposed_fine = new MultiNTensorField3D<T>(numVariables,
        //                                           managementFine,
        //                                           defaultMultiBlockPolicy3D().getBlockCommunicator(),
        //                                           defaultMultiBlockPolicy3D().getCombinedStatistics(),
        //                                           defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>());

        decomposed_t0->periodicity().toggle(0, xPeriodic_);
        decomposed_t12->periodicity().toggle(0, xPeriodic_);
        decomposed_t1->periodicity().toggle(0, xPeriodic_);
        decomposed_fine->periodicity().toggle(0, xPeriodic_);

        decomposed_t0->periodicity().toggle(1, yPeriodic_);
        decomposed_t12->periodicity().toggle(1, yPeriodic_);
        decomposed_t1->periodicity().toggle(1, yPeriodic_);
        decomposed_fine->periodicity().toggle(1, yPeriodic_);

        decomposed_t0->periodicity().toggle(2, zPeriodic_);
        decomposed_t12->periodicity().toggle(2, zPeriodic_);
        decomposed_t1->periodicity().toggle(2, zPeriodic_);
        decomposed_fine->periodicity().toggle(2, zPeriodic_);
    } else {
        decomposed_t0 = 0;
        decomposed_t12 = 0;
        decomposed_t1 = 0;
        decomposed_fine = 0;
    }

    decomposedContainer.push_back(decomposed_t0);
    decomposedContainer.push_back(decomposed_t1);
    decomposedContainer.push_back(decomposed_t12);

    lattice->initialize();
}

template <typename T, template <typename U> class Descriptor>
GridLevelContainer3D<T, Descriptor>::GridLevelContainer3D(
    GridLevelContainer3D<T, Descriptor> const &rhs) :
    dyn(rhs.dyn),
    order(rhs.order),
    level(rhs.level),
    numVariables(rhs.numVariables),
    lattice(rhs.lattice),
    decomposed_t0(rhs.decomposed_t0),
    decomposed_t12(rhs.decomposed_t12),
    decomposed_t1(rhs.decomposed_t1),
    decomposed_fine(rhs.decomposed_fine),
    decomposedContainer(rhs.decomposedContainer),
    parity(rhs.parity),
    owned(rhs.owned)
{ }

template <typename T, template <typename U> class Descriptor>
GridLevelContainer3D<T, Descriptor> &GridLevelContainer3D<T, Descriptor>::operator=(
    GridLevelContainer3D<T, Descriptor> const &rhs)
{
    GridLevelContainer3D<T, Descriptor>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
void GridLevelContainer3D<T, Descriptor>::swap(GridLevelContainer3D<T, Descriptor> &rhs)
{
    std::swap(dyn, rhs.dyn);
    std::swap(order, rhs.order);
    std::swap(level, rhs.level);
    std::swap(numVariables, rhs.numVariables);
    std::swap(lattice, rhs.lattice);
    std::swap(decomposed_t0, rhs.decomposed_t0);
    std::swap(decomposed_t12, rhs.decomposed_t12);
    std::swap(decomposed_t1, rhs.decomposed_t1);
    std::swap(decomposed_fine, rhs.decomposed_fine);
    std::swap(decomposedContainer, rhs.decomposedContainer);
    std::swap(parity, rhs.parity);
    std::swap(owned, rhs.owned);
}

template <typename T, template <typename U> class Descriptor>
GridLevelContainer3D<T, Descriptor>::~GridLevelContainer3D()
{
    if (owned) {
        delete dyn;
        delete lattice;
        delete dataProcessors;
        delete decomposed_t0;
        delete decomposed_t12;
        delete decomposed_t1;
        delete decomposed_fine;
    }
}

// ======================================================================= //
// ====================MultiLevelCoupling3D=============================== //
// ======================================================================= //

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
MultiLevelCoupling3D<T, Descriptor, Engine>::MultiLevelCoupling3D(
    OctreeGridStructure &ogs_, Dynamics<T, Descriptor> *dyn_, plint order_, plint overlapWidth_,
    bool filterAll_, bool xPeriodic_, bool yPeriodic_, bool zPeriodic_) :
    ogs(ogs_), dyn(dyn_), order(order_), overlapWidth(overlapWidth_), filterAll(filterAll_)
{
    // pcout << "Generating grid levels and coupling interfaces." << std::endl;
    generateLevelsAndCouplingInterfaces(xPeriodic_, yPeriodic_, zPeriodic_);
    // pcout << "Initializing tensor fields and integrating processing functionals." << std::endl;
    initializeTensorFields();
    integrateProcessingFunctionals();
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
MultiLevelCoupling3D<T, Descriptor, Engine>::MultiLevelCoupling3D(
    MultiLevelCoupling3D<T, Descriptor, Engine> const &rhs) :
    ogs(rhs.ogs),
    dyn(rhs.dyn),
    order(rhs.order),
    overlapWidth(rhs.overlapWidth),
    gridLevels(rhs.gridLevels),
    couplings(rhs.couplings),
    filterAll(rhs.filterAll)
{ }

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
MultiLevelCoupling3D<T, Descriptor, Engine> &MultiLevelCoupling3D<T, Descriptor, Engine>::operator=(
    MultiLevelCoupling3D<T, Descriptor, Engine> const &rhs)
{
    MultiLevelCoupling3D<T, Descriptor, Engine>(rhs).swap(*this);
    return *this;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelCoupling3D<T, Descriptor, Engine>::swap(
    MultiLevelCoupling3D<T, Descriptor, Engine> &rhs)
{
    std::swap(ogs, rhs.ogs);
    std::swap(dyn, rhs.dyn);
    std::swap(order, rhs.order);
    std::swap(overlapWidth, rhs.overlapWidth);
    std::swap(gridLevels, rhs.gridLevels);
    std::swap(couplings, rhs.couplings);
    std::swap(filterAll, rhs.filterAll);
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
MultiLevelCoupling3D<T, Descriptor, Engine>::~MultiLevelCoupling3D()
{
    delete dyn;
    for (pluint iA = 0; iA < gridLevels.size(); ++iA) {
        delete gridLevels[iA];
    }
    for (pluint iA = 0; iA < couplings.size(); ++iA) {
        delete couplings[iA];
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelCoupling3D<T, Descriptor, Engine>::generateLevelsAndCouplingInterfaces(
    bool xPeriodic, bool yPeriodic, bool zPeriodic)
{
    // computation of the relaxation time vector used for the dynamics
    std::vector<Array<T, Descriptor<T>::q> > relaxationFrequencies(getNumLevels());
    Array<T, Descriptor<T>::q> relFreq = dyn->getRelaxationFrequencies();
    Engine<T, Descriptor> engine;

    for (plint level = 0; level < getNumLevels(); ++level) {
        relaxationFrequencies[level] = relFreq;
        // going from coarse to fine => xDt = 0.5, xInvDx=0 => is useless and
        Array<T, Descriptor<T>::q> relFreqTmp = engine.computeRescaledRelFreq(relFreq, (T)0.5);
        relFreq = relFreqTmp;
    }

    // creation of the coupling interfaces vector
    for (plint level = 0; level < getNumLevels() - 1; ++level) {
        couplings.push_back(new CouplingInterfaces3D(ogs, level + 1, overlapWidth));
    }

    // creation of the lattice and tensor fields
    for (plint level = 0; level < getNumLevels(); ++level) {
        dyn->setRelaxationFrequencies(relaxationFrequencies[level]);
        gridLevels.push_back(new GridLevelContainer3D<T, Descriptor>(
            dyn->clone(), ogs, order, level, xPeriodic, yPeriodic, zPeriodic));
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelCoupling3D<T, Descriptor, Engine>::initialize()
{
    // initialization of the lattice in each grid level
    for (pluint iL = 0; iL < gridLevels.size(); ++iL) {
        gridLevels[iL]->getLattice().initialize();
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelCoupling3D<T, Descriptor, Engine>::initializeTensorFields()
{
    // initialization of the t0 NTenforfields
    for (plint iL = 0; iL < (plint)gridLevels.size() - 1; ++iL) {
        std::vector<boxLogic::DirectedPlane> cToFext_cu = couplings[iL]->getCtoFext_cu();
        for (pluint iA = 0; iA < cToFext_cu.size(); ++iA) {
            applyProcessingFunctional(
                new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>(
                    (T)0.5, T(), gridLevels[iL]->getOrder()),
                cToFext_cu[iA].bb, gridLevels[iL]->getLattice(), gridLevels[iL]->getDecomposedT0());
            applyProcessingFunctional(
                new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>(
                    (T)0.5, T(), gridLevels[iL]->getOrder()),
                cToFext_cu[iA].bb, gridLevels[iL]->getLattice(), gridLevels[iL]->getDecomposedT1());
        }
        std::vector<Box3D> fToCint_fu = couplings[iL]->getFtoCint_fu();
        for (pluint iA = 0; iA < fToCint_fu.size(); ++iA) {
            applyProcessingFunctional(
                new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>(
                    (T)2, T(), gridLevels[iL + 1]->getOrder()),
                fToCint_fu[iA], gridLevels[iL + 1]->getLattice(),
                gridLevels[iL]->getDecomposedFine());
        }
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelCoupling3D<T, Descriptor, Engine>::integrateProcessingFunctionals()
{
    // initialization of the t0 NTenforfields
    for (pluint iL = 0; iL < gridLevels.size(); ++iL) {
        integrateProcessingFunctionals(iL);
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelCoupling3D<T, Descriptor, Engine>::integrateCopyAndSpatialInterpolation(
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
void MultiLevelCoupling3D<T, Descriptor, Engine>::integrateProcessingFunctionals(plint iL)
{
    if (iL < (plint)(gridLevels.size() - 1)) {
        // 1. Decompose cTensor at t1
        std::vector<boxLogic::DirectedPlane> cToFext = couplings[iL]->getCtoFext_cu();
        for (plint iA = 0; iA < (plint)cToFext.size(); ++iA) {
            integrateProcessingFunctional(
                new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>(
                    (T)0.5, T(), gridLevels[iL]->getOrder()),
                cToFext[iA].bb, gridLevels[iL]->getLattice(), gridLevels[iL]->getDecomposedT1(),
                lattice_decomposeT1);
            integrateProcessingFunctional(
                new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>(
                    (T)0.5, T(), gridLevels[iL]->getOrder()),
                cToFext[iA].bb, gridLevels[iL]->getLattice(), gridLevels[iL]->getDecomposedT0(),
                lattice_decomposeT0);
        }

        for (pluint iA = 0; iA < cToFext.size(); ++iA) {
            integrateProcessingFunctional(
                new TemporalInterpolationFunctional3D<T>(), cToFext[iA].bb,
                gridLevels[iL]->getDecomposedContainer(), container_timeInterpolation);
            std::vector<MultiBlock3D *> swappedArg(gridLevels[iL]->getDecomposedContainer());
            std::swap(swappedArg[0], swappedArg[1]);
            integrateProcessingFunctional(
                new TemporalInterpolationFunctional3D<T>(), cToFext[iA].bb, swappedArg,
                container_timeInterpolation);
        }

        integrateCopyAndSpatialInterpolation(
            gridLevels[iL]->getDecomposedT12(), gridLevels[iL]->getDecomposedFine(), *couplings[iL],
            decomposedT12_spatialInterpolation);

        std::vector<boxLogic::Plane> cToFint = couplings[iL]->getCtoFint_fu();
        for (pluint iA = 0; iA < cToFint.size(); ++iA) {
            integrateProcessingFunctional(
                new RecomposeFunctional3D<T, Descriptor>(gridLevels[iL + 1]->getOrder()),
                cToFint[iA].bb, gridLevels[iL + 1]->getLattice(),
                gridLevels[iL]->getDecomposedFine(), lattice_recomposeFine);
        }

        integrateCopyAndSpatialInterpolation(
            gridLevels[iL]->getDecomposedT0(), gridLevels[iL]->getDecomposedFine(), *couplings[iL],
            decomposedT1_spatialInterpolation);
        integrateCopyAndSpatialInterpolation(
            gridLevels[iL]->getDecomposedT1(), gridLevels[iL]->getDecomposedFine(), *couplings[iL],
            decomposedT1_spatialInterpolation);

        // fine to coarse coupling
        std::vector<Box3D> fToCint_fu = couplings[iL]->getFtoCint_fu();
        for (pluint iA = 0; iA < fToCint_fu.size(); ++iA) {
            integrateProcessingFunctional(
                new DecomposeAndRescaleFunctional3D<T, Descriptor, Engine>(
                    (T)2, T(), gridLevels[iL + 1]->getOrder()),
                fToCint_fu[iA], gridLevels[iL + 1]->getLattice(),
                gridLevels[iL]->getDecomposedFine(), lattice_decomposeFine);
        }

        std::vector<boxLogic::Plane> fToCint_cu = couplings[iL]->getFtoCint_cu();
        for (pluint iA = 0; iA < fToCint_cu.size(); ++iA) {
            integrateProcessingFunctional(
                new CopyAndFilterNonEquilibriumFunctional3D<T, Descriptor>(filterAll),
                fToCint_cu[iA].bb, gridLevels[iL]->getDecomposedT12(),
                gridLevels[iL]->getDecomposedFine(), decomposedT12_copyAndFilter);
        }

        std::vector<boxLogic::Edge> fToCbcEdges_cu = couplings[iL]->getFtoCboundaryEdges_cu();
        for (pluint iA = 0; iA < fToCbcEdges_cu.size(); ++iA) {
            integrateProcessingFunctional(
                new CopyFunctional3D<T, Descriptor>(), fToCbcEdges_cu[iA].bb,
                gridLevels[iL]->getDecomposedT12(), gridLevels[iL]->getDecomposedFine(),
                decomposedT12_copyBc);
        }

        for (pluint iA = 0; iA < fToCint_cu.size(); ++iA) {
            integrateProcessingFunctional(
                new RecomposeFunctional3D<T, Descriptor>(gridLevels[iL]->getOrder()),
                fToCint_cu[iA].bb, gridLevels[iL]->getLattice(), gridLevels[iL]->getDecomposedT12(),
                lattice_recomposeCoarse);
        }

        for (pluint iA = 0; iA < fToCbcEdges_cu.size(); ++iA) {
            integrateProcessingFunctional(
                new RecomposeFunctional3D<T, Descriptor>(gridLevels[iL]->getOrder()),
                fToCbcEdges_cu[iA].bb, gridLevels[iL]->getLattice(),
                gridLevels[iL]->getDecomposedT12(), lattice_recomposeBc);
        }
    }
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void MultiLevelCoupling3D<T, Descriptor, Engine>::collideAndStream(
    plint iL, std::map<plint, bool> const &useExecuteInternalProcessors,
    std::vector<plint> const &extProcIds, bool computeStats, plint statsId,
    std::vector<std::vector<plint> > const &ids,
    std::vector<std::vector<std::vector<T> > > &results)
{
    std::map<plint, bool>::const_iterator it = useExecuteInternalProcessors.find(iL);
    PLB_ASSERT(it != useExecuteInternalProcessors.end());
    global::timer("gr_collideAndStream").start();
    global::timer("gr_collideAndStream_" + util::val2str(iL)).start();
    if (it->second) {
        gridLevels[iL]->getLattice().executeInternalProcessors();  // collision coarse
        gridLevels[iL]->getLattice().incrementTime();  // Increment the time manually, because,
                                                       // although it is done implicitly by
                                                       // collideAndStream(), it cannot be invoked
                                                       // from a data processor.
    } else {
        gridLevels[iL]->getLattice().collideAndStream();  // collision coarse
    }
    global::timer("gr_collideAndStream_" + util::val2str(iL)).stop();
    global::timer("gr_collideAndStream").stop();

    if (iL < (plint)(gridLevels.size() - 1)) {
        global::timer("gr_coupling").start();
        global::timer("gr_coupling_" + util::val2str(iL)).start();
        // 1. Decompose cTensor at t1
        global::timer("gr_decomposeAndRescale_" + util::val2str(iL)).start();
        global::timer("gr_decomposeAndRescale").start();
        if (gridLevels[iL]->getParity() == 0) {
            gridLevels[iL]->getLattice().executeInternalProcessors(lattice_decomposeT1);
        } else {
            gridLevels[iL]->getLattice().executeInternalProcessors(lattice_decomposeT0);
        }
        global::timer("gr_decomposeAndRescale").stop();
        global::timer("gr_decomposeAndRescale_" + util::val2str(iL)).stop();
        // 2. Interpolate cTensor at t12
        global::timer("gr_temporalInterpolate").start();
        global::timer("gr_temporalInterpolate_" + util::val2str(iL)).start();
        gridLevels[iL]->getDecomposedContainer()[0]->executeInternalProcessors(
            container_timeInterpolation, false);
        global::timer("gr_temporalInterpolate").stop();
        global::timer("gr_temporalInterpolate_" + util::val2str(iL)).stop();

        global::timer("gr_coupling_" + util::val2str(iL)).stop();
        global::timer("gr_coupling").stop();

        collideAndStream(
            iL + 1, useExecuteInternalProcessors, extProcIds, computeStats, statsId, ids,
            results);  // collision fine
        global::timer("gr_coupling_" + util::val2str(iL)).start();
        global::timer("gr_coupling").start();

        // Coarse to fine coupling:
        global::timer("gr_copyAndSpatialInterpolation").start();
        global::timer("gr_copyAndSpatialInterpolation_" + util::val2str(iL)).start();
        gridLevels[iL]->getDecomposedT12().executeInternalProcessors(
            decomposedT12_spatialInterpolation);
        global::timer("gr_copyAndSpatialInterpolation").stop();
        global::timer("gr_copyAndSpatialInterpolation_" + util::val2str(iL)).stop();

        global::timer("gr_recompose_fine").start();
        gridLevels[iL + 1]->getLattice().executeInternalProcessors(lattice_recomposeFine, false);
        global::timer("gr_recompose_fine").stop();
        for (plint iA = 0; iA < (plint)extProcIds.size(); ++iA) {
            gridLevels[iL + 1]->getLattice().executeInternalProcessors(extProcIds[iA], false);
        }

        if (computeStats) {
            gridLevels[iL + 1]->getLattice().executeInternalProcessors(statsId, false);
            for (pluint iA = 0; iA < ids[iL + 1].size(); ++iA) {
                results[iL + 1][iA].push_back(
                    gridLevels[iL + 1]->getLattice().getInternalStatistics().getSum(
                        ids[iL + 1][iA]));
            }
        }

        global::timer("gr_coupling").stop();
        global::timer("gr_coupling_" + util::val2str(iL)).stop();
        collideAndStream(
            iL + 1, useExecuteInternalProcessors, extProcIds, computeStats, statsId, ids,
            results);  // collision fine
        global::timer("gr_coupling_" + util::val2str(iL)).start();
        global::timer("gr_coupling").start();

        global::timer("gr_copyAndSpatialInterpolation_" + util::val2str(iL)).start();
        global::timer("gr_copyAndSpatialInterpolation").start();
        gridLevels[iL]->getDecomposedT1().executeInternalProcessors(
            decomposedT1_spatialInterpolation);
        global::timer("gr_copyAndSpatialInterpolation").stop();
        global::timer("gr_copyAndSpatialInterpolation_" + util::val2str(iL)).stop();

        global::timer("gr_recompose_fine").start();
        gridLevels[iL + 1]->getLattice().executeInternalProcessors(lattice_recomposeFine, false);
        global::timer("gr_recompose_fine").stop();
        for (plint iA = 0; iA < (plint)extProcIds.size(); ++iA) {
            gridLevels[iL + 1]->getLattice().executeInternalProcessors(extProcIds[iA]);
        }
        if (computeStats) {
            gridLevels[iL + 1]->getLattice().executeInternalProcessors(statsId, false);
            for (pluint iA = 0; iA < ids[iL + 1].size(); ++iA) {
                results[iL + 1][iA].push_back(
                    gridLevels[iL + 1]->getLattice().getInternalStatistics().getSum(
                        ids[iL + 1][iA]));
            }
        }

        // fine to coarse coupling
        std::vector<Box3D> fToCint_fu = couplings[iL]->getFtoCint_fu();
        global::timer("gr_decomposeAndRescale").start();
        global::timer("gr_decomposeAndRescale_" + util::val2str(iL)).start();
        gridLevels[iL + 1]->getLattice().executeInternalProcessors(lattice_decomposeFine);
        global::timer("gr_decomposeAndRescale_" + util::val2str(iL)).stop();
        global::timer("gr_decomposeAndRescale").stop();

        global::timer("gr_copyAndFilter").start();
        gridLevels[iL]->getDecomposedT12().executeInternalProcessors(decomposedT12_copyAndFilter);
        global::timer("gr_copyAndFilter").stop();

        global::timer("gr_copy").start();
        gridLevels[iL]->getDecomposedT12().executeInternalProcessors(decomposedT12_copyBc);
        global::timer("gr_copy").stop();

        global::timer("gr_recompose_coarse").start();

        gridLevels[iL]->getLattice().executeInternalProcessors(lattice_recomposeCoarse, false);
        global::timer("gr_recompose_coarse").stop();
        global::timer("gr_recompose_bc").start();
        gridLevels[iL]->getLattice().executeInternalProcessors(lattice_recomposeBc, false);
        global::timer("gr_recompose_bc").stop();

        if (iL == 0) {
            for (plint iA = 0; iA < (plint)extProcIds.size(); ++iA) {
                gridLevels[iL]->getLattice().executeInternalProcessors(extProcIds[iA], false);
            }
            if (computeStats) {
                gridLevels[iL]->getLattice().executeInternalProcessors(statsId, false);
                for (pluint iA = 0; iA < ids[iL].size(); ++iA) {
                    results[iL][iA].push_back(
                        gridLevels[iL]->getLattice().getInternalStatistics().getSum(ids[iL][iA]));
                }
            }
        }

        // swaps decomposedT0 and decomposedT1
        gridLevels[iL]->swapTimeSteps();
        global::timer("gr_coupling").stop();
        global::timer("gr_coupling_" + util::val2str(iL)).stop();
    }
}

/*
template<typename T, template <typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
void MultiLevelCoupling3D<T,Descriptor,Engine>::collideAndStream(plint iL, std::map<plint, bool>
const& useExecuteInternalProcessors)
{
    std::map<plint, bool>::const_iterator it = useExecuteInternalProcessors.find(iL);
    PLB_ASSERT(it != useExecuteInternalProcessors.end());
    global::timer("gr_collideAndStream_"+util::val2str(iL)).start();
    if (it->second) {
        gridLevels[iL]->getLattice().executeInternalProcessors();  // collision coarse
    } else {
        gridLevels[iL]->getLattice().collideAndStream();  // collision coarse
    }
    global::timer("gr_collideAndStream_"+util::val2str(iL)).stop();

    if (iL < (plint)(gridLevels.size()-1)) {
        // 1. Decompose cTensor at t1
        global::timer("gr_decomposeAndRescale").start();
        gridLevels[iL]->getLattice().executeInternalProcessors(decomposeCoarseT1);
        global::timer("gr_decomposeAndRescale").stop();

        // 2. Interpolate cTensor at t12
        global::timer("gr_temporalInterpolate").start();
        gridLevels[iL]->getLattice().executeInternalProcessors(timeInterpolation);
        global::timer("gr_temporalInterpolate").stop();

        collideAndStream(iL+1, useExecuteInternalProcessors);  // collision fine
        // Coarse to fine coupling:
        global::timer("gr_copyAndSpatialInterpolation").start();
        gridLevels[iL]->getLattice().executeInternalProcessors(spatialInterpolationT12);
        global::timer("gr_copyAndSpatialInterpolation").stop();

        global::timer("gr_recompose").start();
        gridLevels[iL+1]->getLattice().executeInternalProcessors(recomposeFineT12);
        global::timer("gr_recompose").stop();

        collideAndStream(iL+1, useExecuteInternalProcessors);  // collision fine

        global::timer("gr_copyAndSpatialInterpolation").start();
        gridLevels[iL]->getLattice().executeInternalProcessors(spatialInterpolationT1);
        global::timer("gr_copyAndSpatialInterpolation").stop();

        global::timer("gr_recompose").start();
        gridLevels[iL+1]->getLattice().executeInternalProcessors(recomposeFineT1);
        global::timer("gr_recompose").stop();

        // fine to coarse coupling
        global::timer("gr_decomposeAndRescale").start();
        gridLevels[iL+1]->getLattice().executeInternalProcessors(decomposeFine);
        global::timer("gr_decomposeAndRescale").stop();

        global::timer("gr_copyAndFilter").start();
        gridLevels[iL]->getLattice().executeInternalProcessors(copyAndFilter);
        global::timer("gr_copyAndFilter").stop();

        global::timer("gr_copy").start();
        gridLevels[iL]->getLattice().executeInternalProcessors(copyBc);
        global::timer("gr_copy").stop();

        global::timer("gr_recompose").start();
        gridLevels[iL]->getLattice().executeInternalProcessors(recomposeCoarse);
        gridLevels[iL]->getLattice().executeInternalProcessors(recomposeBc);
        global::timer("gr_recompose").stop();

        // swaps decomposedT0 and decomposedT1
        gridLevels[iL]->swapTimeSteps();
    }
}
*/

}  // namespace plb

#endif  // COUPLING_INTERFACE_GENERATOR_3D_HH
