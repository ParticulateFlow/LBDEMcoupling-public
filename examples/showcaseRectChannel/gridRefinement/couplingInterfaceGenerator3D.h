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

#ifndef COUPLING_INTERFACE_GENERATOR_3D_H
#define COUPLING_INTERFACE_GENERATOR_3D_H

#include <map>
#include <vector>

#include "core/globalDefs.h"
#include "gridRefinement/boxLogic3D.h"
#include "gridRefinement/multiLevel3D.h"
#include "gridRefinement/octreeGridStructure.h"
#include "multiBlock/group3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

class CouplingInterfaces3D {
public:
    CouplingInterfaces3D(OctreeGridStructure const &ogs_, int level_, int overlapWidth);

    void writeInterfaces(double dx, const Array<double, 3> &pos) const;

    // coarse to fine interfaces
    std::vector<boxLogic::DirectedPlane> const &getCtoFext_cu() const
    {
        return coarseToFinePlanesCoarseUnitsExtended;
    }
    std::vector<Box3D> boxCtoFext_cu() const
    {
        std::vector<Box3D> result(coarseToFinePlanesCoarseUnitsExtended.size());
        for (pluint i = 0; i < coarseToFinePlanesCoarseUnitsExtended.size(); ++i) {
            result[i] = coarseToFinePlanesCoarseUnitsExtended[i].bb;
        }
        return result;
    }
    std::vector<boxLogic::DirectedPlane> const &getCtoFplanes_cu() const
    {
        return coarseToFinePlanesCoarseUnits;
    }
    std::vector<boxLogic::DirectedEdge> const &getCtoFbulkEdges_cu() const
    {
        return coarseToFineEdgesCoarseUnits;
    }
    std::vector<boxLogic::DirectedCorner> const &getCtoFbulkCorners_cu() const
    {
        return coarseToFineCornersCoarseUnits;
    }
    std::vector<boxLogic::DirectedEdge> const &getCtoFboundaryEdges_cu() const
    {
        return coarseToFineBoundaryEdgesCoarseUnits;
    }
    std::vector<boxLogic::DirectedCorner> const &getCtoFboundaryCorners_cu() const
    {
        return coarseToFineBoundaryCornersCoarseUnits;
    }
    std::vector<boxLogic::Plane> const &getCtoFint_fu() const
    {
        return coarseToFinePlanesFineUnits;
    }
    std::vector<Box3D> boxCtoFint_fu() const
    {
        std::vector<Box3D> result(coarseToFinePlanesFineUnits.size());
        for (pluint i = 0; i < coarseToFinePlanesFineUnits.size(); ++i) {
            result[i] = coarseToFinePlanesFineUnits[i].bb;
        }
        return result;
    }

    // fine to coarse interfaces
    std::vector<Box3D> const &getFtoCint_fu() const
    {
        return fineToCoarsePlanesFineUnits;
    }
    std::vector<boxLogic::Plane> const &getFtoCint_cu() const
    {
        return fineToCoarsePlanesCoarseUnits;
    }
    std::vector<Box3D> boxFtoCint_cu() const
    {
        std::vector<Box3D> result(fineToCoarsePlanesCoarseUnits.size());
        for (pluint i = 0; i < fineToCoarsePlanesCoarseUnits.size(); ++i) {
            result[i] = fineToCoarsePlanesCoarseUnits[i].bb;
        }
        return result;
    }
    std::vector<boxLogic::Edge> const &getFtoCboundaryEdges_cu() const
    {
        return fineToCoarseBoundaryEdgesCoarseUnits;
    }
    std::vector<Box3D> boxFtoCboundaryEdges_cu() const
    {
        std::vector<Box3D> result(fineToCoarseBoundaryEdgesCoarseUnits.size());
        for (pluint i = 0; i < fineToCoarseBoundaryEdgesCoarseUnits.size(); ++i) {
            result[i] = fineToCoarseBoundaryEdgesCoarseUnits[i].bb;
        }
        return result;
    }

private:
    OctreeGridStructure const &ogs;
    int level;

private:
    std::vector<boxLogic::Plane>
        coarseToFinePlanesFineUnits;  // used for recompose of fine interface
    std::vector<boxLogic::DirectedPlane>
        coarseToFinePlanesCoarseUnits;  // used for interpolation Processing Functionals (PF)
    std::vector<boxLogic::DirectedEdge> coarseToFineEdgesCoarseUnits;  // used for interpolation PF
    std::vector<boxLogic::DirectedCorner>
        coarseToFineCornersCoarseUnits;  // used for interpolation PF
    std::vector<boxLogic::DirectedPlane>
        coarseToFinePlanesCoarseUnitsExtended;       // used for decompose of coarse interface
    std::vector<Box3D> fineToCoarsePlanesFineUnits;  // used for decompose of fine interfare
    std::vector<boxLogic::Plane>
        fineToCoarsePlanesCoarseUnits;  // used for filtering and recompose of coarse interface

    std::vector<boxLogic::DirectedEdge>
        coarseToFineBoundaryEdgesCoarseUnits;  // special inepolation for BC edge
    std::vector<boxLogic::DirectedCorner>
        coarseToFineBoundaryCornersCoarseUnits;  // special inepolation for BC corner
    std::vector<boxLogic::Edge> fineToCoarseBoundaryEdgesCoarseUnits;  // used for interpolation PF
};

// This class conatins a grid level. A locally refined grid usually contains many of these
// grid levels.
template <typename T, template <typename U> class Descriptor>
class GridLevelContainer3D {
public:
    GridLevelContainer3D(
        Dynamics<T, Descriptor> *dyn_, const OctreeGridStructure &ogs, plint order_, plint level_,
        bool xPeriodic_, bool yPeriodic_, bool zPeriodic_);
    GridLevelContainer3D(std::vector<Group3D *> &groups, plint order_, plint level_);

    ~GridLevelContainer3D();

    MultiContainerBlock3D &getDataProcessors()
    {
        return *dataProcessors;
    }
    MultiBlockLattice3D<T, Descriptor> &getLattice() const
    {
        return *lattice;
    }
    MultiNTensorField3D<T> &getDecomposedT0() const
    {
        return *decomposed_t0;
    }
    MultiNTensorField3D<T> &getDecomposedT12() const
    {
        return *decomposed_t12;
    }
    MultiNTensorField3D<T> &getDecomposedT1() const
    {
        return *decomposed_t1;
    }
    MultiNTensorField3D<T> &getDecomposedFine() const
    {
        return *decomposed_fine;
    }
    std::vector<MultiBlock3D *> &getDecomposedContainer()
    {
        return decomposedContainer;
    }

    plint getParity() const
    {
        return parity;
    }

    void swapTimeSteps()
    {
        std::swap(decomposed_t0, decomposed_t1);
        std::swap(decomposedContainer[0], decomposedContainer[1]);
        parity = 1 - parity;
    }

    plint const &getOrder() const
    {
        return order;
    }

private:
    GridLevelContainer3D(GridLevelContainer3D<T, Descriptor> const &rhs);
    GridLevelContainer3D<T, Descriptor> &operator=(GridLevelContainer3D<T, Descriptor> const &rhs);
    void swap(GridLevelContainer3D<T, Descriptor> &rhs);

private:
    Dynamics<T, Descriptor> *dyn;
    plint order, level, numVariables;

private:
    MultiContainerBlock3D *dataProcessors;
    MultiBlockLattice3D<T, Descriptor> *lattice;
    MultiNTensorField3D<T> *decomposed_t0, *decomposed_t12, *decomposed_t1, *decomposed_fine;
    std::vector<MultiBlock3D *> decomposedContainer;
    plint parity;
    bool owned;
};

// Contains all the necessary information for a complete grid refined lattice.
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
class MultiLevelCoupling3D : public MultiLevel3D {
public:
    // The dynamics dyn must be right for the coarsest instantiated level.
    MultiLevelCoupling3D(
        OctreeGridStructure &ogs_, Dynamics<T, Descriptor> *dyn_, plint order_,
        plint overlapWidth_ = 1, bool filterAll_ = true, bool xPeriodic_ = false,
        bool yPeriodic_ = false, bool zPeriodic_ = false);
    MultiLevelCoupling3D(MultiLevelCoupling3D<T, Descriptor, Engine> const &rhs);
    MultiLevelCoupling3D<T, Descriptor, Engine> &operator=(
        MultiLevelCoupling3D<T, Descriptor, Engine> const &rhs);
    void swap(MultiLevelCoupling3D<T, Descriptor, Engine> &rhs);
    ~MultiLevelCoupling3D();

    void initialize();

    void writeInterfaces(double dx, const Array<double, 3> &pos) const
    {
        for (plint iL = 0; iL < (plint)couplings.size(); ++iL) {
            couplings[iL]->writeInterfaces(dx, pos);
        }
    }

    void collideAndStream(
        plint iL, std::map<plint, bool> const &useExecuteInternalProcessors,
        std::vector<plint> const &extProcIds, bool computeStats, plint statsId,
        std::vector<std::vector<plint> > const &ids,
        std::vector<std::vector<std::vector<T> > > &results);

    virtual MultiBlockLattice3D<T, Descriptor> const &getLevel(plint iL) const
    {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return gridLevels[iL]->getLattice();
    }

    virtual MultiBlockLattice3D<T, Descriptor> &getLevel(plint iL)
    {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return gridLevels[iL]->getLattice();
    }

    GridLevelContainer3D<T, Descriptor> const &getGridLevelContainer(plint iL) const
    {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return *gridLevels[iL];
    }

    GridLevelContainer3D<T, Descriptor> &getGridLevelContainer(plint iL)
    {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return *gridLevels[iL];
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

private:
    void generateLevelsAndCouplingInterfaces(bool xPeriodic, bool yPeriodic, bool zPeriodic);
    void integrateProcessingFunctionals();
    void integrateProcessingFunctionals(plint iL);
    void integrateCopyAndSpatialInterpolation(
        MultiNTensorField3D<T> &cTensor, MultiNTensorField3D<T> &fTensor,
        CouplingInterfaces3D &coupling, plint numExec);

private:
    OctreeGridStructure ogs;
    Dynamics<T, Descriptor> *dyn;
    plint order;
    plint overlapWidth;

private:
    std::vector<GridLevelContainer3D<T, Descriptor> *> gridLevels;
    std::vector<CouplingInterfaces3D *> couplings;
    bool filterAll;
};

}  // namespace plb

#endif  // COUPLING_INTERFACE_GENERATOR_3D_H
