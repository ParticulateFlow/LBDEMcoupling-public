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
 * Various factories that use a multiGridManagement2D -- Header file
 */

#ifndef MULTI_GRID_GENERATOR_2D_HH
#define MULTI_GRID_GENERATOR_2D_HH

#include "multiGrid/multiGridGenerator2D.h"

namespace plb {

/// Use the MultiGridManagement2D to generate a verctor of lattices that represent the multi grid.
template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlockLattice2D<T, Descriptor> *> generateLattices(
    MultiGridManagement2D management, std::vector<Dynamics<T, Descriptor> *> backgroundDynamics,
    std::vector<BlockCommunicator2D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics, plint envelopeWidth)
{
    // get the MPI id's (if available)
    std::vector<std::vector<plint> > mpiProcesses = management.getMpiProcesses();
    std::vector<std::vector<Box2D> > bulks = management.getBulks();
    PLB_PRECONDITION(mpiProcesses.size() == 0 || mpiProcesses.size() == bulks.size());

    // allocation of the result vector
    std::vector<MultiBlockLattice2D<T, Descriptor> *> multiBlocks(bulks.size());
    for (pluint iLevel = 0; iLevel < bulks.size(); ++iLevel) {
        // create a thread attribution
        ExplicitThreadAttribution *threadAttribution = new ExplicitThreadAttribution;
        PLB_ASSERT(mpiProcesses.size() == 0 || mpiProcesses[iLevel].size() == bulks[iLevel].size());
        SparseBlockStructure2D geometry(management.getBoundingBox(iLevel));
        for (pluint iBlock = 0; iBlock < bulks[iLevel].size(); ++iBlock) {
            plint mpiProcess = mpiProcesses.size() > 0 ? mpiProcesses[iLevel][iBlock] : 0;
            plint blockId = geometry.nextIncrementalId();
            Box2D bulk = bulks[iLevel][iBlock];
            geometry.addBlock(bulk, bulk, blockId);
            threadAttribution->addBlock(blockId, mpiProcess);
        }
        MultiBlockManagement2D blockManagement(geometry, threadAttribution, envelopeWidth, iLevel);
        multiBlocks[iLevel] = new MultiBlockLattice2D<T, Descriptor>(
            blockManagement, communicators[iLevel], combinedStatistics[iLevel],
            defaultMultiBlockPolicy2D().getMultiCellAccess<T, Descriptor>(),
            backgroundDynamics[iLevel]);
        multiBlocks[iLevel]->periodicity().toggleAll(false);
    }
    return multiBlocks;
}

template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlockLattice2D<T, Descriptor> *> generateLattices(
    MultiGridManagement2D management, std::vector<Dynamics<T, Descriptor> *> backgroundDynamics,
    plint envelopeWidth)
{
    return generateLattices(
        management, backgroundDynamics,
        defaultMultiGridPolicy2D().getBlockCommunicator<T>(management.getNumLevels()),
        defaultMultiGridPolicy2D().getCombinedStatistics(management.getNumLevels()), envelopeWidth);
}

/** Given a MultiGridManagement2D object, create the MultiScalarField2D that
 *   represent the hierarchy.
 */
template <typename T>
std::vector<MultiScalarField2D<T> *> generateScalarFields(
    MultiGridManagement2D const &management, std::vector<BlockCommunicator2D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics)
{
    std::vector<std::vector<Box2D> > bulks = management.getBulks();
    std::vector<std::vector<plint> > mpiProcesses = management.getMpiProcesses();

    PLB_PRECONDITION(mpiProcesses.size() == 0 || mpiProcesses.size() == bulks.size());
    plint envelopeWidth = 1;
    std::vector<MultiScalarField2D<T> *> scalars(bulks.size());
    for (plint iLevel = 0; iLevel < (plint)bulks.size(); ++iLevel) {
        PLB_ASSERT(mpiProcesses.size() == 0 || mpiProcesses[iLevel].size() == bulks[iLevel].size());
        ExplicitThreadAttribution *threadAttribution = new ExplicitThreadAttribution;
        SparseBlockStructure2D geometry(management.getBoundingBox(iLevel));
        for (pluint iBlock = 0; iBlock < bulks[iLevel].size(); ++iBlock) {
            plint mpiProcess = mpiProcesses.size() > 0 ? mpiProcesses[iLevel][iBlock] : 0;
            plint blockId = geometry.nextIncrementalId();
            geometry.addBlock(bulks[iLevel][iBlock], bulks[iLevel][iBlock], blockId);
            threadAttribution->addBlock(blockId, mpiProcess);
        }
        MultiBlockManagement2D blockManagement(geometry, threadAttribution, envelopeWidth, iLevel);
        scalars[iLevel] = new MultiScalarField2D<T>(
            blockManagement, communicators[iLevel], combinedStatistics[iLevel],
            defaultMultiBlockPolicy2D().getMultiScalarAccess<T>());
    }

    return scalars;
}

template <typename T, int nDim>
std::vector<MultiTensorField2D<T, nDim> *> generateTensorFields(
    MultiGridManagement2D const &management, std::vector<BlockCommunicator2D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics)
{
    std::vector<MultiTensorField2D<T, nDim> *> fields(management.getNumLevels());
    std::vector<std::vector<Box2D> > bulks = management.getBulks();
    std::vector<std::vector<plint> > ids = management.getMpiProcesses();
    PLB_PRECONDITION(ids.size() == 0 || ids.size() == bulks.size());

    plint envelopeWidth = 1;

    fields.resize(bulks.size());
    for (plint iLevel = 0; iLevel < (plint)bulks.size(); ++iLevel) {
        ExplicitThreadAttribution *threadAttribution = new ExplicitThreadAttribution;
        PLB_ASSERT(ids.size() == 0 || ids[iLevel].size() == bulks[iLevel].size());
        SparseBlockStructure2D geometry(management.getBoundingBox(iLevel));
        for (pluint iBlock = 0; iBlock < bulks[iLevel].size(); ++iBlock) {
            plint mpiProcess = ids.size() > 0 ? ids[iLevel][iBlock] : 0;
            plint blockId = geometry.nextIncrementalId();
            geometry.addBlock(bulks[iLevel][iBlock], bulks[iLevel][iBlock], blockId);
            threadAttribution->addBlock(blockId, mpiProcess);
        }

        MultiBlockManagement2D blockManagement(geometry, threadAttribution, envelopeWidth, iLevel);
        fields[iLevel] = new MultiTensorField2D<T, nDim>(
            blockManagement, communicators[iLevel], combinedStatistics[iLevel],
            defaultMultiBlockPolicy2D().getMultiTensorAccess<T, nDim>());
    }

    return fields;
}

template <typename T, template <typename U> class Descriptor>
void LinearInterpolationFineGridInterfaceInstantiator<T, Descriptor>::instantiateDataProcessors(
    Box2D fineGridInterface, MultiBlockLattice2D<T, Descriptor> &coarseLattice,
    MultiBlockLattice2D<T, Descriptor> &fineLattice, plint direction, plint orientation)

{
    PLB_PRECONDITION(fineGridInterface.getNx() == 1 || fineGridInterface.getNy() == 1);
    plint scalingOrder = 0;  // decompose in rho, u and fneq

    RescaleEngine<T, Descriptor> *rescaleEngine =
        new ConvectiveRescaleEngine<T, Descriptor>(scalingOrder);

    // computation of the direction according to the refinement
    plint GRdirection = fineGridInterface.getNy() == 1 ? 0 : 1;
    // computation of the direction according to Palabos BC convention
    plint BCdirection = direction;

    Box2D reducedFineGridInterface(fineGridInterface);
    if (GRdirection == 0) {  // Interface extends in x-direction.
        reducedFineGridInterface.x0 += 1;
        reducedFineGridInterface.x1 -= 1;
    } else {  // interface extends in y-direction.
        reducedFineGridInterface.y0 += 1;
        reducedFineGridInterface.y1 -= 1;
    }

    Box2D lowerFineEnd, upperFineEnd;
    lowerFineEnd.x0 = lowerFineEnd.x1 = fineGridInterface.x0;
    lowerFineEnd.y0 = lowerFineEnd.y1 = fineGridInterface.y0;
    upperFineEnd.x0 = upperFineEnd.x1 = fineGridInterface.x1;
    upperFineEnd.y0 = upperFineEnd.y1 = fineGridInterface.y1;

    // Number of time steps executed by the fine grid during a coarse iteration.
    plint numTimeSteps = 2;
    plint executionTime = numTimeSteps - 1;

    setCompositeDynamics(
        fineLattice, fineGridInterface.multiply(2),
        new FineGridBoundaryDynamics<T, Descriptor>(
            new NoDynamics<T, Descriptor>, fineLattice.getTimeCounter(), numTimeSteps,
            rescaleEngine->getDecompositionOrder()));

    // The coarse processors are at processing level 1, because they need to access information
    //   from within the coarse envelope.
    plint processorLevelCoarse = 1;

    // Copy, rescale, and interpolate from coarse to fine in the bulk of the interface.
    integrateProcessingFunctional(
        new CopyCoarseToFineLinearInterp2D<T, Descriptor, Descriptor>(
            rescaleEngine->clone(), BCdirection, orientation),
        reducedFineGridInterface, coarseLattice, fineLattice, processorLevelCoarse);

    // Copy, rescale, and interpolate from coarse to fine at the lower end of the interface.
    integrateProcessingFunctional(
        new CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor, Descriptor>(
            rescaleEngine->clone(), BCdirection, orientation, 1),
        lowerFineEnd, coarseLattice, fineLattice, processorLevelCoarse);

    // Copy, rescale, and interpolate from coarse to fine at the upper end of the interface.
    integrateProcessingFunctional(
        new CopyCoarseToFineBoundaryLinearInterp2D<T, Descriptor, Descriptor>(
            rescaleEngine->clone(), BCdirection, orientation, -1),
        upperFineEnd, coarseLattice, fineLattice, processorLevelCoarse);

    // Add a data processor which imposes a time-cyclic behavior on the fine grid
    //   boundary-dynamics.
    integrateProcessingFunctional(
        new Copy_t1_to_t0_2D<T, Descriptor>(numTimeSteps, executionTime),
        fineGridInterface.multiply(2), fineLattice);

    delete rescaleEngine;
}

template <typename T, template <typename U> class Descriptor>
void CubicInterpolationFineGridInterfaceInstantiator<T, Descriptor>::instantiateDataProcessors(
    Box2D fineGridInterface, MultiBlockLattice2D<T, Descriptor> &coarseLattice,
    MultiBlockLattice2D<T, Descriptor> &fineLattice, plint direction, plint orientation)

{
    PLB_PRECONDITION(fineGridInterface.getNx() == 1 || fineGridInterface.getNy() == 1);
    plint scalingOrder = 0;  // to decompose in rho,u and fneq
    RescaleEngine<T, Descriptor> *rescaleEngine =
        new ConvectiveRescaleEngine<T, Descriptor>(scalingOrder);

    // computation of the direction according to the refinement
    plint GRdirection = fineGridInterface.getNy() == 1 ? 0 : 1;
    // computation of the direction according to Palabos BC convention
    plint BCdirection = direction;

    Box2D reducedFineGridInterface(fineGridInterface);
    if (GRdirection == 0) {  // Interface extends in x-direction.
        reducedFineGridInterface.x0 += 1;
        reducedFineGridInterface.x1 -= 1;
    } else {  // interface extends in y-direction.
        reducedFineGridInterface.y0 += 1;
        reducedFineGridInterface.y1 -= 1;
    }

    Box2D lowerFineEnd, upperFineEnd;
    lowerFineEnd.x0 = lowerFineEnd.x1 = fineGridInterface.x0;
    lowerFineEnd.y0 = lowerFineEnd.y1 = fineGridInterface.y0;
    upperFineEnd.x0 = upperFineEnd.x1 = fineGridInterface.x1;
    upperFineEnd.y0 = upperFineEnd.y1 = fineGridInterface.y1;

    // Number of time steps executed by the fine grid during a coarse iteration.
    plint numTimeSteps = 2;
    plint executionTime = numTimeSteps - 1;

    setCompositeDynamics(
        fineLattice, fineGridInterface.multiply(2),
        new FineGridBoundaryDynamics<T, Descriptor>(
            new NoDynamics<T, Descriptor>, fineLattice.getTimeCounter(), numTimeSteps,
            rescaleEngine->getDecompositionOrder()));

    // The coarse processors are at processing level 1, because they need to access information
    //   from within the coarse envelope.
    plint processorLevelCoarse = 1;

    // Copy, rescale, and interpolate from coarse to fine in the bulk of the interface.
    integrateProcessingFunctional(
        new CopyCoarseToFineCubicInterp2D<T, Descriptor, Descriptor>(
            rescaleEngine->clone(), BCdirection, orientation),
        reducedFineGridInterface, coarseLattice, fineLattice, processorLevelCoarse);

    // define where to interpolate for the boundaries
    plint toRight = 1;
    plint toLeft = -1;

    // Copy, rescale, and interpolate from coarse to fine at the lower end of the interface.
    integrateProcessingFunctional(
        new CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor, Descriptor>(
            rescaleEngine->clone(), BCdirection, orientation, toRight),
        lowerFineEnd, coarseLattice, fineLattice, processorLevelCoarse);

    // Copy, rescale, and interpolate from coarse to fine at the upper end of the interface.
    integrateProcessingFunctional(
        new CopyCoarseToFineBoundaryCubicInterp2D<T, Descriptor, Descriptor>(
            rescaleEngine->clone(), BCdirection, orientation, toLeft),
        upperFineEnd, coarseLattice, fineLattice, processorLevelCoarse);

    // Add helpers to complete the interpolations in sites near the corners
    Box2D reducedLowerFineGridBorder(
        reducedFineGridInterface.x0, reducedFineGridInterface.x0, reducedFineGridInterface.y0,
        reducedFineGridInterface.y0);
    Box2D reducedUpperFineGridBorder(
        reducedFineGridInterface.x1, reducedFineGridInterface.x1, reducedFineGridInterface.y1,
        reducedFineGridInterface.y1);

    integrateProcessingFunctional(
        new CopyCoarseToFineBoundaryHelper2D<T, Descriptor, Descriptor>(
            rescaleEngine->clone(), BCdirection, orientation, toRight),
        reducedLowerFineGridBorder, coarseLattice, fineLattice, processorLevelCoarse);

    integrateProcessingFunctional(
        new CopyCoarseToFineBoundaryHelper2D<T, Descriptor, Descriptor>(
            rescaleEngine->clone(), BCdirection, orientation, toLeft),
        reducedUpperFineGridBorder, coarseLattice, fineLattice, processorLevelCoarse);

    // Add a data processor which imposes a time-cyclic behavior on the fine grid
    //   boundary-dynamics.
    integrateProcessingFunctional(
        new Copy_t1_to_t0_2D<T, Descriptor>(numTimeSteps, executionTime),
        fineGridInterface.multiply(2), fineLattice);

    delete rescaleEngine;
}

template <typename T, template <typename U> class Descriptor>
void FilteredCoarseGridInterfaceInstantiator<T, Descriptor>::instantiateDataProcessors(
    Box2D coarseGridInterface, MultiBlockLattice2D<T, Descriptor> &coarseLattice,
    MultiBlockLattice2D<T, Descriptor> &fineLattice, plint direction, plint orientation)
{
    // Number of time steps executed by the fine grid during a coarse iteration.
    plint numTimeSteps = 2;
    plint executionTime = numTimeSteps - 1;
    // Add a data processor which copies and rescales from the fine to the coarse lattice.
    plint scalingOrder = 0;  // to decompose in rho,u and fneq

    integrateProcessingFunctional(
        new CopyFineToCoarseWithFiltering2D<T, Descriptor, Descriptor>(
            new ConvectiveRescaleEngine<T, Descriptor>(scalingOrder), numTimeSteps, executionTime,
            direction, orientation),
        coarseGridInterface.multiply(2), fineLattice, coarseLattice, 2);
}

template <typename T, template <typename U> class Descriptor>
void UnfilteredCoarseGridInterfaceInstantiator<T, Descriptor>::instantiateDataProcessors(
    Box2D coarseGridInterface, MultiBlockLattice2D<T, Descriptor> &coarseLattice,
    MultiBlockLattice2D<T, Descriptor> &fineLattice, plint direction, plint orientation)
{
    PLB_PRECONDITION(coarseGridInterface.getNx() == 1 || coarseGridInterface.getNy() == 1);

    // Number of time steps executed by the fine grid during a coarse iteration.
    plint numTimeSteps = 2;
    plint executionTime = numTimeSteps - 1;
    // Add a data processor which copies and rescales from the fine to the coarse lattice.
    plint scalingOrder = 0;  // to decompose in rho,u and fneq

    integrateProcessingFunctional(
        new CopyFineToCoarse2D<T, Descriptor, Descriptor>(
            new ConvectiveRescaleEngine<T, Descriptor>(scalingOrder), numTimeSteps, executionTime,
            direction, orientation),
        coarseGridInterface.multiply(2), fineLattice, coarseLattice, 2);
}

template <typename T, template <typename U> class Descriptor>

std::unique_ptr<MultiGridLattice2D<T, Descriptor> >
    MultiGridGenerator2D<T, Descriptor>::createRefinedLatticeCubicInterpolationNoFiltering(
        MultiGridManagement2D management, Dynamics<T, Descriptor> *backgroundDynamics,
        plint behaviorLevel)
{
    return std::unique_ptr<MultiGridLattice2D<T, Descriptor> >(
        new MultiGridLattice2D<T, Descriptor>(
            management,
            defaultMultiGridPolicy2D().getBlockCommunicator<T>(management.getNumLevels()),
            defaultMultiGridPolicy2D().getCombinedStatistics(management.getNumLevels()),
            backgroundDynamics, behaviorLevel,
            new CubicInterpolationFineGridInterfaceInstantiator<T, Descriptor>(),
            new UnfilteredCoarseGridInterfaceInstantiator<T, Descriptor>()));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridLattice2D<T, Descriptor> >
    MultiGridGenerator2D<T, Descriptor>::createRefinedLatticeLinearInterpolationNoFiltering(
        MultiGridManagement2D management, Dynamics<T, Descriptor> *backgroundDynamics,
        plint behaviorLevel)
{
    return std::unique_ptr<MultiGridLattice2D<T, Descriptor> >(
        new MultiGridLattice2D<T, Descriptor>(
            management,
            defaultMultiGridPolicy2D().getBlockCommunicator<T>(management.getNumLevels()),
            defaultMultiGridPolicy2D().getCombinedStatistics(management.getNumLevels()),
            backgroundDynamics, behaviorLevel,
            new LinearInterpolationFineGridInterfaceInstantiator<T, Descriptor>(),
            new UnfilteredCoarseGridInterfaceInstantiator<T, Descriptor>()));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridLattice2D<T, Descriptor> >
    MultiGridGenerator2D<T, Descriptor>::createRefinedLatticeCubicInterpolationFiltering(
        MultiGridManagement2D management, Dynamics<T, Descriptor> *backgroundDynamics,
        plint behaviorLevel)
{
    return std::unique_ptr<MultiGridLattice2D<T, Descriptor> >(
        new MultiGridLattice2D<T, Descriptor>(
            management,
            defaultMultiGridPolicy2D().getBlockCommunicator<T>(management.getNumLevels()),
            defaultMultiGridPolicy2D().getCombinedStatistics(management.getNumLevels()),
            backgroundDynamics, behaviorLevel,
            new CubicInterpolationFineGridInterfaceInstantiator<T, Descriptor>(),
            new FilteredCoarseGridInterfaceInstantiator<T, Descriptor>()));
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiGridLattice2D<T, Descriptor> >
    MultiGridGenerator2D<T, Descriptor>::createRefinedLatticeLinearInterpolationFiltering(
        MultiGridManagement2D management, Dynamics<T, Descriptor> *backgroundDynamics,
        plint behaviorLevel)
{
    return std::unique_ptr<MultiGridLattice2D<T, Descriptor> >(
        new MultiGridLattice2D<T, Descriptor>(
            management,
            defaultMultiGridPolicy2D().getBlockCommunicator<T>(management.getNumLevels()),
            defaultMultiGridPolicy2D().getCombinedStatistics(management.getNumLevels()),
            backgroundDynamics, behaviorLevel,
            new LinearInterpolationFineGridInterfaceInstantiator<T, Descriptor>(),
            new FilteredCoarseGridInterfaceInstantiator<T, Descriptor>()));
}

}  // namespace plb

#endif  // MULTI_GRID_GENERATOR_2D_HH
