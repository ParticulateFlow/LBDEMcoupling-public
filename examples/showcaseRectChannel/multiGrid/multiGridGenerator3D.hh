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
 * Various factories that use a multiGridManagement3D -- Header file
 */

#ifndef MULTI_GRID_GENERATOR_3D_HH
#define MULTI_GRID_GENERATOR_3D_HH

#include "multiGrid/coarseGridProcessors3D.h"
#include "multiGrid/domainDivision3D.h"
#include "multiGrid/fineGridProcessors3D.h"
#include "multiGrid/helperFineGridProcessors3D.h"
#include "multiGrid/multiGridGenerator3D.h"

namespace plb {

Box3D computeCopyReducedBulk(Box3D domain)
{
    Box3D result(domain);
    if (domain.getNx() == 1) {
        result.y0++;
        result.y1--;
        result.z0++;
        result.z1--;
    }
    if (domain.getNy() == 1) {
        result.x0++;
        result.x1--;
        result.z0++;
        result.z1--;
    }
    if (domain.getNz() == 1) {
        result.x0++;
        result.x1--;
        result.y0++;
        result.y1--;
    }
    return result;
}

void computeCopyEdges(Box3D domain, std::vector<Box3D> &edges)
{
    edges.resize(4);
    if (domain.getNx() == 1) {
        edges[0] = Box3D(domain.x0, domain.x1, domain.y0 + 1, domain.y1 - 1, domain.z0, domain.z0);
        edges[1] = Box3D(domain.x0, domain.x1, domain.y0 + 1, domain.y1 - 1, domain.z1, domain.z1);
        edges[2] = Box3D(domain.x0, domain.x1, domain.y0, domain.y0, domain.z0 + 1, domain.z1 - 1);
        edges[3] = Box3D(domain.x0, domain.x1, domain.y1, domain.y1, domain.z0 + 1, domain.z1 - 1);
    }
    if (domain.getNy() == 1) {
        edges[0] = Box3D(domain.x0 + 1, domain.x1 - 1, domain.y0, domain.y1, domain.z0, domain.z0);
        edges[1] = Box3D(domain.x0 + 1, domain.x1 - 1, domain.y0, domain.y1, domain.z1, domain.z1);
        edges[2] = Box3D(domain.x0, domain.x0, domain.y0, domain.y1, domain.z0 + 1, domain.z1 - 1);
        edges[3] = Box3D(domain.x1, domain.x1, domain.y0, domain.y1, domain.z0 + 1, domain.z1 - 1);
    }
    if (domain.getNz() == 1) {
        edges[0] = Box3D(domain.x0 + 1, domain.x1 - 1, domain.y0, domain.y0, domain.z0, domain.z1);
        edges[1] = Box3D(domain.x0 + 1, domain.x1 - 1, domain.y1, domain.y1, domain.z0, domain.z1);
        edges[2] = Box3D(domain.x0, domain.x0, domain.y0 + 1, domain.y1 - 1, domain.z0, domain.z1);
        edges[3] = Box3D(domain.x1, domain.x1, domain.y0 + 1, domain.y1 - 1, domain.z0, domain.z1);
    }
}

template <typename T>
void computeFilteringIndicesEdges(Box3D domain, std::vector<std::vector<plint> > &indices)
{
    indices.resize(4);
    if (domain.getNx() == 1) {
        std::vector<plint> index1 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 2, 1>();
        //         std::vector<plint> index2 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 2, 0>();
        indices[0].insert(indices[0].end(), index1.begin(), index1.end());
        //         indices[0].insert(indices[0].end(),index2.begin(),index2.end());

        std::vector<plint> index3 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 2, -1>();
        //         std::vector<plint> index4 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 2,  0>();
        indices[1].insert(indices[1].end(), index3.begin(), index3.end());
        //         indices[1].insert(indices[1].end(),index4.begin(),index4.end());

        std::vector<plint> index5 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 1, 1>();
        //         std::vector<plint> index6 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 1, 0>();
        indices[2].insert(indices[2].end(), index5.begin(), index5.end());
        //         indices[2].insert(indices[2].end(),index6.begin(),index6.end());

        std::vector<plint> index7 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 1, -1>();
        //         std::vector<plint> index8 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 1,  0>();
        indices[3].insert(indices[3].end(), index7.begin(), index7.end());
        //         indices[3].insert(indices[3].end(),index8.begin(),index8.end());
    }

    if (domain.getNy() == 1) {
        std::vector<plint> index1 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 2, 1>();
        //         std::vector<plint> index2 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 0, 0>();
        indices[0].insert(indices[0].end(), index1.begin(), index1.end());
        //         indices[0].insert(indices[0].end(),index2.begin(),index2.end());

        std::vector<plint> index3 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 2, -1>();
        //         std::vector<plint> index4 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 0, 0>();
        indices[1].insert(indices[1].end(), index3.begin(), index3.end());
        //         indices[1].insert(indices[1].end(),index4.begin(),index4.end());

        std::vector<plint> index5 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 0, 1>();
        //         std::vector<plint> index6 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 2, 0>();
        indices[2].insert(indices[2].end(), index5.begin(), index5.end());
        //         indices[2].insert(indices[2].end(),index6.begin(),index6.end());

        std::vector<plint> index7 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 0, -1>();
        //         std::vector<plint> index8 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 2, 0>();
        indices[3].insert(indices[3].end(), index7.begin(), index7.end());
        //         indices[3].insert(indices[3].end(),index8.begin(),index8.end());
    }

    if (domain.getNz() == 1) {
        std::vector<plint> index1 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 1, 1>();
        //         std::vector<plint> index2 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 0, 0>();
        indices[0].insert(indices[0].end(), index1.begin(), index1.end());
        //         indices[0].insert(indices[0].end(),index2.begin(),index2.end());

        std::vector<plint> index3 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 1, -1>();
        //         std::vector<plint> index4 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 0,  0>();
        indices[1].insert(indices[1].end(), index3.begin(), index3.end());
        //         indices[1].insert(indices[1].end(),index4.begin(),index4.end());

        std::vector<plint> index5 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 0, 1>();
        //         std::vector<plint> index6 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 1, 0>();
        indices[2].insert(indices[2].end(), index5.begin(), index5.end());
        //         indices[2].insert(indices[2].end(),index6.begin(),index6.end());

        std::vector<plint> index7 =
            indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 0, -1>();
        //         std::vector<plint> index8 =
        //         indexTemplates::subIndex<descriptors::D3Q27Descriptor<T>, 1,  0>();
        indices[3].insert(indices[3].end(), index7.begin(), index7.end());
        //         indices[3].insert(indices[3].end(),index8.begin(),index8.end());
    }
}

void computeCopyCorners(Box3D domain, std::vector<Box3D> &corners)
{
    corners.resize(4);
    if (domain.getNx() == 1) {
        corners[0] = Box3D(domain.x0, domain.x1, domain.y0, domain.y0, domain.z0, domain.z0);
        corners[1] = Box3D(domain.x0, domain.x1, domain.y0, domain.y0, domain.z1, domain.z1);
        corners[2] = Box3D(domain.x0, domain.x1, domain.y1, domain.y1, domain.z0, domain.z0);
        corners[3] = Box3D(domain.x0, domain.x1, domain.y1, domain.y1, domain.z1, domain.z1);
    }
    if (domain.getNy() == 1) {
        corners[0] = Box3D(domain.x0, domain.x0, domain.y0, domain.y1, domain.z0, domain.z0);
        corners[1] = Box3D(domain.x0, domain.x0, domain.y0, domain.y1, domain.z1, domain.z1);
        corners[2] = Box3D(domain.x1, domain.x1, domain.y0, domain.y1, domain.z0, domain.z0);
        corners[3] = Box3D(domain.x1, domain.x1, domain.y0, domain.y1, domain.z1, domain.z1);
    }
    if (domain.getNz() == 1) {
        corners[0] = Box3D(domain.x0, domain.x0, domain.y0, domain.y0, domain.z0, domain.z1);
        corners[1] = Box3D(domain.x0, domain.x0, domain.y1, domain.y1, domain.z0, domain.z1);
        corners[2] = Box3D(domain.x1, domain.x1, domain.y0, domain.y0, domain.z0, domain.z1);
        corners[3] = Box3D(domain.x1, domain.x1, domain.y1, domain.y1, domain.z0, domain.z1);
    }
}

/// Use the MultiGridManagement3D to generate a vector of lattices that represent the multi grid.
template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlockLattice3D<T, Descriptor> *> generateLattices(
    MultiGridManagement3D management, std::vector<Dynamics<T, Descriptor> *> backgroundDynamics,
    std::vector<BlockCommunicator3D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics)
{
    // get the MPI id's (if available)
    std::vector<std::vector<plint> > mpiProcesses = management.getMpiProcesses();
    std::vector<std::vector<Box3D> > bulks = management.getBulks();
    PLB_PRECONDITION(mpiProcesses.size() == 0 || mpiProcesses.size() == bulks.size());
    // for cubic interpolation
    plint envelopeWidth = 2;
    // allocation of the result vector
    std::vector<MultiBlockLattice3D<T, Descriptor> *> multiBlocks(bulks.size());
    for (pluint iLevel = 0; iLevel < bulks.size(); ++iLevel) {
        // create a thread attribution
        ExplicitThreadAttribution *threadAttribution = new ExplicitThreadAttribution;
        PLB_ASSERT(mpiProcesses.size() == 0 || mpiProcesses[iLevel].size() == bulks[iLevel].size());
        SparseBlockStructure3D geometry(management.getBoundingBox(iLevel));
        for (pluint iBlock = 0; iBlock < bulks[iLevel].size(); ++iBlock) {
            plint mpiProcess = mpiProcesses.size() > 0 ? mpiProcesses[iLevel][iBlock] : 0;
            plint blockId = geometry.nextIncrementalId();
            Box3D bulk = bulks[iLevel][iBlock];
            geometry.addBlock(bulk, bulk, blockId);  // TODO: add uniqueBulk
            threadAttribution->addBlock(blockId, mpiProcess);
        }
        MultiBlockManagement3D blockManagement(geometry, threadAttribution, envelopeWidth, iLevel);
        multiBlocks[iLevel] = new MultiBlockLattice3D<T, Descriptor>(
            blockManagement, communicators[iLevel], combinedStatistics[iLevel],
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
            backgroundDynamics[iLevel]);
        multiBlocks[iLevel]->periodicity().toggleAll(false);
    }
    return multiBlocks;
}

template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlockLattice3D<T, Descriptor> *> generateLattices(
    MultiGridManagement3D management, std::vector<Dynamics<T, Descriptor> *> backgroundDynamics)
{
    return generateLattices(
        management, backgroundDynamics,
        defaultMultiGridPolicy3D().getBlockCommunicator(management.getNumLevels()),
        defaultMultiGridPolicy3D().getCombinedStatistics(management.getNumLevels()));
}

template <typename T, template <typename U> class Descriptor>
void createInterfaces(
    std::vector<MultiBlockLattice3D<T, Descriptor> *> &multiBlocks,
    MultiGridManagement3D management)
{
    std::vector<std::vector<Box3D> > fineGridInterfaces = management.getFineInterface();
    std::vector<std::vector<Box3D> > coarseGridInterfaces = management.getCoarseInterface();

    // Create the part of the coupling in which a lattice acts as a coarse lattice.
    for (pluint iCoarse = 0; iCoarse < multiBlocks.size() - 1; ++iCoarse) {
        for (pluint iInterface = 0; iInterface < coarseGridInterfaces[iCoarse].size(); ++iInterface)
        {
            createCoarseGridInterface(
                iCoarse, coarseGridInterfaces[iCoarse][iInterface], multiBlocks);
        }
    }
    // Create the part of the coupling in which a lattice acts as a fine lattice.
    for (pluint iFine = 1; iFine < multiBlocks.size(); ++iFine) {
        for (pluint iInterface = 0; iInterface < fineGridInterfaces[iFine].size(); ++iInterface) {
            createFineGridInterface(iFine - 1, fineGridInterfaces[iFine][iInterface], multiBlocks);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void createCoarseGridInterface(
    plint coarseLevel, Box3D coarseGridInterface,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> &multiBlocks)
{
    PLB_PRECONDITION(
        coarseGridInterface.getNx() == 1 || coarseGridInterface.getNy() == 1
        || coarseGridInterface.getNz() == 1);

    pluint fineLevel = coarseLevel + 1;
    PLB_PRECONDITION(fineLevel < multiBlocks.size());

    MultiBlockLattice3D<T, Descriptor> &coarseLattice = *multiBlocks[coarseLevel];
    MultiBlockLattice3D<T, Descriptor> &fineLattice = *multiBlocks[fineLevel];

    // Number of time steps executed by the fine grid during a coarse iteration.
    plint numTimeSteps = 2;
    plint executionTime = numTimeSteps - 1;
    // Add a data processor which copies and rescales from the fine to the coarse lattice.
    plint scalingOrder = 0;

    /// ****** NO FILTERING **********/
    /*    integrateProcessingFunctional (
                new CopyFineToCoarse3D<T, Descriptor,Descriptor> (
                    new ConvectiveRescaleEngine<T,Descriptor>(scalingOrder), numTimeSteps,
       executionTime), coarseGridInterface.multiply(2), fineLattice, coarseLattice );
    */

    /// ******* WITH EDGES AND CORNERS EXPLICITLY ********/
    /*     Box3D reducedBulk = computeCopyReducedBulk(coarseGridInterface);
         std::vector<Box3D> edges;
         computeCopyEdges(coarseGridInterface, edges);

         std::vector<Box3D> corners;
         computeCopyCorners(coarseGridInterface, corners);

         std::vector<plint> bulkIndices;
         std::vector<std::vector<plint> > edgesFilteringIndices(edges.size());
         std::vector<std::vector<plint> > cornersFilteringIndices(corners.size());


         computeFilteringIndicesEdges<T>(coarseGridInterface, edgesFilteringIndices);

         // always start from 1, cause the 0 is already and index
         for (plint iPop=1; iPop<descriptors::D3Q27Descriptor<T>::q; ++iPop){
             bulkIndices.push_back(iPop);
         }

        // assign the copy corresponding to the reduced bulk (full indices)
         integrateProcessingFunctional (
                 new CopyFineToCoarseWithFiltering3D<T, Descriptor,Descriptor> (
                     new ConvectiveRescaleEngine<T,Descriptor>(scalingOrder),
                     numTimeSteps, executionTime, bulkIndices ),
                 reducedBulk.multiply(2),
                 fineLattice, coarseLattice, 2 ); // I (Helen) added the 2 here...

        // assign the copy corresponding to the edges (reduced indices)
         for (pluint iEdge=0; iEdge<edges.size(); ++iEdge){
             integrateProcessingFunctional (
                 new CopyFineToCoarseWithFiltering3D<T, Descriptor,Descriptor> (
                     new ConvectiveRescaleEngine<T,Descriptor>(scalingOrder),
                     numTimeSteps, executionTime, edgesFilteringIndices[iEdge] ),
                 edges[iEdge].multiply(2),
                 fineLattice, coarseLattice, 2 ); // I (Helen) added the 2 here...
         }

        // finally do the same for the corners
         for (pluint iCorner=0; iCorner<corners.size(); ++iCorner){
             integrateProcessingFunctional (
                 new CopyFineToCoarseWithFiltering3D<T, Descriptor,Descriptor> (
                     new ConvectiveRescaleEngine<T,Descriptor>(scalingOrder),
                     numTimeSteps, executionTime, cornersFilteringIndices[iCorner] ),
                 corners[iCorner].multiply(2),
                 fineLattice, coarseLattice, 2 ); // I (Helen) added the 2 here...
         }
    */
    /// ***** ALL AT ONCE - ORIGINAL ********/
    /*     std::vector<plint> bulkIndices;
         // always start from 1, cause the 0 is already and index
         for (plint iPop=1; iPop<descriptors::D3Q27Descriptor<T>::q; ++iPop){
             bulkIndices.push_back(iPop);
         }

         integrateProcessingFunctional (
                 new CopyFineToCoarseWithFiltering3D<T, Descriptor,Descriptor> (
                     new ConvectiveRescaleEngine<T,Descriptor>(scalingOrder),
                     numTimeSteps, executionTime, bulkIndices ),
                 coarseGridInterface.multiply(2),
                 fineLattice, coarseLattice );
    */
    /// ****** ALL AT ONCE - MODIFIED *******/
    std::vector<plint> bulkIndices;
    // always start from 1, cause the 0 is already an index
    for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
        bulkIndices.push_back(iPop);
    }

    integrateProcessingFunctional(
        new CopyFineToCoarseWithFiltering3D<T, Descriptor, Descriptor>(
            new ConvectiveRescaleEngine<T, Descriptor>(scalingOrder), numTimeSteps, executionTime,
            bulkIndices),
        coarseGridInterface.multiply(2), fineLattice, coarseLattice,
        2);  // I (Helen) added the 2 here...
}

template <typename T, template <typename U> class Descriptor>
void createFineGridInterface(
    plint coarseLevel, Box3D fineGridInterface,
    std::vector<MultiBlockLattice3D<T, Descriptor> *> &multiBlocks)
{
    plint scalingOrder = 0;
    RescaleEngine<T, Descriptor> *rescaleEngine =
        new ConvectiveRescaleEngine<T, Descriptor>(scalingOrder);

    PLB_PRECONDITION(
        fineGridInterface.getNx() == 1 || fineGridInterface.getNy() == 1
        || fineGridInterface.getNz() == 1);
    pluint fineLevel = coarseLevel + 1;
    PLB_PRECONDITION(fineLevel < multiBlocks.size());

    MultiBlockLattice3D<T, Descriptor> &coarseLattice = *multiBlocks[coarseLevel];
    MultiBlockLattice3D<T, Descriptor> &fineLattice = *multiBlocks[fineLevel];

    plint direction;
    if (fineGridInterface.getNx() == 1)
        direction = 0;
    if (fineGridInterface.getNy() == 1)
        direction = 1;
    if (fineGridInterface.getNz() == 1)
        direction = 2;

    // Number of time steps executed by the fine grid during a coarse iteration.
    plint numTimeSteps = 2;
    plint executionTime = numTimeSteps - 1;

    // Instantiate specific time-interpolation dynamics on the interface of the fine grid.
    setCompositeDynamics(
        fineLattice, fineGridInterface.multiply(2),
        new FineGridBoundaryDynamics<T, Descriptor>(
            new NoDynamics<T, Descriptor>, fineLattice.getTimeCounter(), numTimeSteps,
            rescaleEngine->getDecompositionOrder()));

    // The coarse processors are at processing level 1, because they need to access information
    //   from within the coarse envelope.
    plint processorLevelCoarse = 1;

    // Add interpolation processors
    // find out the orientation of the plane
    if (direction == 0) {  // Plan YZ
        Box3D reducedCoarseInterface(
            fineGridInterface.x0, fineGridInterface.x1, fineGridInterface.y0 + 1,
            fineGridInterface.y1 - 1, fineGridInterface.z0 + 1, fineGridInterface.z1 - 1);

        // rewritten data processor for the bulk interpolation over the fine grid
        integrateProcessingFunctional(
            new ScalarCubicInterpolationYZ<T, Descriptor>(rescaleEngine->clone()),
            reducedCoarseInterface, coarseLattice, fineLattice, processorLevelCoarse);

        // rewritten data processor for the edge interpolation over the fine grid
        Box3D edge1(
            fineGridInterface.x0, fineGridInterface.x1, fineGridInterface.y0 + 1,
            fineGridInterface.y1 - 1, fineGridInterface.z0, fineGridInterface.z0);
        Box3D edge2(
            fineGridInterface.x0, fineGridInterface.x1, fineGridInterface.y0 + 1,
            fineGridInterface.y1 - 1, fineGridInterface.z1, fineGridInterface.z1);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationYZLineY3D<T, Descriptor>(1, rescaleEngine->clone()), edge1,
            coarseLattice, fineLattice, processorLevelCoarse);
        integrateProcessingFunctional(
            new ScalarCubicInterpolationYZLineY3D<T, Descriptor>(-1, rescaleEngine->clone()), edge2,
            coarseLattice, fineLattice, processorLevelCoarse);

        Box3D edge3(
            fineGridInterface.x0, fineGridInterface.x1, fineGridInterface.y0, fineGridInterface.y0,
            fineGridInterface.z0 + 1, fineGridInterface.z1 - 1);
        Box3D edge4(
            fineGridInterface.x0, fineGridInterface.x1, fineGridInterface.y1, fineGridInterface.y1,
            fineGridInterface.z0 + 1, fineGridInterface.z1 - 1);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationYZLineZ3D<T, Descriptor>(1, rescaleEngine->clone()), edge3,
            coarseLattice, fineLattice, processorLevelCoarse);
        integrateProcessingFunctional(
            new ScalarCubicInterpolationYZLineZ3D<T, Descriptor>(-1, rescaleEngine->clone()), edge4,
            coarseLattice, fineLattice, processorLevelCoarse);

        // rewritten data processor for the corner interpolation over the fine grid
        Box3D cornerNN(
            fineGridInterface.x0, fineGridInterface.x1, fineGridInterface.y0, fineGridInterface.y0,
            fineGridInterface.z0, fineGridInterface.z0);
        Box3D cornerNP(
            fineGridInterface.x0, fineGridInterface.x1, fineGridInterface.y0, fineGridInterface.y0,
            fineGridInterface.z1, fineGridInterface.z1);
        Box3D cornerPN(
            fineGridInterface.x0, fineGridInterface.x1, fineGridInterface.y1, fineGridInterface.y1,
            fineGridInterface.z0, fineGridInterface.z0);
        Box3D cornerPP(
            fineGridInterface.x0, fineGridInterface.x1, fineGridInterface.y1, fineGridInterface.y1,
            fineGridInterface.z1, fineGridInterface.z1);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationYZCorner3D<T, Descriptor>(1, 1, rescaleEngine->clone()),
            cornerNN, coarseLattice, fineLattice, processorLevelCoarse);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationYZCorner3D<T, Descriptor>(1, -1, rescaleEngine->clone()),
            cornerNP, coarseLattice, fineLattice, processorLevelCoarse);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationYZCorner3D<T, Descriptor>(-1, 1, rescaleEngine->clone()),
            cornerPN, coarseLattice, fineLattice, processorLevelCoarse);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationYZCorner3D<T, Descriptor>(-1, -1, rescaleEngine->clone()),
            cornerPP, coarseLattice, fineLattice, processorLevelCoarse);
    }

    if (direction == 1) {  // Plan XZ
        Box3D reducedCoarseInterface(
            fineGridInterface.x0 + 1, fineGridInterface.x1 - 1, fineGridInterface.y0,
            fineGridInterface.y1, fineGridInterface.z0 + 1, fineGridInterface.z1 - 1);

        // rewritten data processor for the bulk interpolation over the fine grid
        integrateProcessingFunctional(
            new ScalarCubicInterpolationXZ<T, Descriptor>(rescaleEngine->clone()),
            reducedCoarseInterface, coarseLattice, fineLattice, processorLevelCoarse);

        // rewritten data processor for the edge interpolation over the fine grid
        Box3D edge1(
            fineGridInterface.x0 + 1, fineGridInterface.x1 - 1, fineGridInterface.y0,
            fineGridInterface.y1, fineGridInterface.z0, fineGridInterface.z0);
        Box3D edge2(
            fineGridInterface.x0 + 1, fineGridInterface.x1 - 1, fineGridInterface.y0,
            fineGridInterface.y1, fineGridInterface.z1, fineGridInterface.z1);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXZLineX3D<T, Descriptor>(1, rescaleEngine->clone()), edge1,
            coarseLattice, fineLattice, processorLevelCoarse);
        integrateProcessingFunctional(
            new ScalarCubicInterpolationXZLineX3D<T, Descriptor>(-1, rescaleEngine->clone()), edge2,
            coarseLattice, fineLattice, processorLevelCoarse);

        Box3D edge3(
            fineGridInterface.x0, fineGridInterface.x0, fineGridInterface.y0, fineGridInterface.y1,
            fineGridInterface.z0 + 1, fineGridInterface.z1 - 1);
        Box3D edge4(
            fineGridInterface.x1, fineGridInterface.x1, fineGridInterface.y0, fineGridInterface.y1,
            fineGridInterface.z0 + 1, fineGridInterface.z1 - 1);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXZLineZ3D<T, Descriptor>(1, rescaleEngine->clone()), edge3,
            coarseLattice, fineLattice, processorLevelCoarse);
        integrateProcessingFunctional(
            new ScalarCubicInterpolationXZLineZ3D<T, Descriptor>(-1, rescaleEngine->clone()), edge4,
            coarseLattice, fineLattice, processorLevelCoarse);

        Box3D cornerNN(
            fineGridInterface.x0, fineGridInterface.x0, fineGridInterface.y0, fineGridInterface.y1,
            fineGridInterface.z0, fineGridInterface.z0);
        Box3D cornerNP(
            fineGridInterface.x0, fineGridInterface.x0, fineGridInterface.y0, fineGridInterface.y1,
            fineGridInterface.z1, fineGridInterface.z1);
        Box3D cornerPN(
            fineGridInterface.x1, fineGridInterface.x1, fineGridInterface.y0, fineGridInterface.y1,
            fineGridInterface.z0, fineGridInterface.z0);
        Box3D cornerPP(
            fineGridInterface.x1, fineGridInterface.x1, fineGridInterface.y0, fineGridInterface.y1,
            fineGridInterface.z1, fineGridInterface.z1);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXZCorner3D<T, Descriptor>(1, 1, rescaleEngine->clone()),
            cornerNN, coarseLattice, fineLattice, processorLevelCoarse);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXZCorner3D<T, Descriptor>(1, -1, rescaleEngine->clone()),
            cornerNP, coarseLattice, fineLattice, processorLevelCoarse);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXZCorner3D<T, Descriptor>(-1, 1, rescaleEngine->clone()),
            cornerPN, coarseLattice, fineLattice, processorLevelCoarse);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXZCorner3D<T, Descriptor>(-1, -1, rescaleEngine->clone()),
            cornerPP, coarseLattice, fineLattice, processorLevelCoarse);
    }

    if (direction == 2) {  // Plan XY
        Box3D reducedCoarseInterface(
            fineGridInterface.x0 + 1, fineGridInterface.x1 - 1, fineGridInterface.y0 + 1,
            fineGridInterface.y1 - 1, fineGridInterface.z0, fineGridInterface.z1);

        // rewritten data processor for the bulk interpolation over the fine grid
        integrateProcessingFunctional(
            new ScalarCubicInterpolationXY<T, Descriptor>(rescaleEngine->clone()),
            reducedCoarseInterface, coarseLattice, fineLattice, processorLevelCoarse);

        // rewritten data processor for the edge interpolation over the fine grid
        Box3D edge1(
            fineGridInterface.x0 + 1, fineGridInterface.x1 - 1, fineGridInterface.y0,
            fineGridInterface.y0, fineGridInterface.z0, fineGridInterface.z1);
        Box3D edge2(
            fineGridInterface.x0 + 1, fineGridInterface.x1 - 1, fineGridInterface.y1,
            fineGridInterface.y1, fineGridInterface.z0, fineGridInterface.z1);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXYLineX3D<T, Descriptor>(1, rescaleEngine->clone()), edge1,
            coarseLattice, fineLattice, processorLevelCoarse);
        integrateProcessingFunctional(
            new ScalarCubicInterpolationXYLineX3D<T, Descriptor>(-1, rescaleEngine->clone()), edge2,
            coarseLattice, fineLattice, processorLevelCoarse);

        Box3D edge3(
            fineGridInterface.x0, fineGridInterface.x0, fineGridInterface.y0 + 1,
            fineGridInterface.y1 - 1, fineGridInterface.z0, fineGridInterface.z1);
        Box3D edge4(
            fineGridInterface.x1, fineGridInterface.x1, fineGridInterface.y0 + 1,
            fineGridInterface.y1 - 1, fineGridInterface.z0, fineGridInterface.z1);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXYLineY3D<T, Descriptor>(1, rescaleEngine->clone()), edge3,
            coarseLattice, fineLattice, processorLevelCoarse);
        integrateProcessingFunctional(
            new ScalarCubicInterpolationXYLineY3D<T, Descriptor>(-1, rescaleEngine->clone()), edge4,
            coarseLattice, fineLattice, processorLevelCoarse);

        // rewritten data processor for the corner interpolation over the fine grid
        Box3D cornerNN(
            fineGridInterface.x0, fineGridInterface.x0, fineGridInterface.y0, fineGridInterface.y0,
            fineGridInterface.z0, fineGridInterface.z1);
        Box3D cornerNP(
            fineGridInterface.x0, fineGridInterface.x0, fineGridInterface.y1, fineGridInterface.y1,
            fineGridInterface.z0, fineGridInterface.z0);
        Box3D cornerPN(
            fineGridInterface.x1, fineGridInterface.x1, fineGridInterface.y0, fineGridInterface.y0,
            fineGridInterface.z0, fineGridInterface.z1);
        Box3D cornerPP(
            fineGridInterface.x1, fineGridInterface.x1, fineGridInterface.y1, fineGridInterface.y1,
            fineGridInterface.z0, fineGridInterface.z1);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXYCorner3D<T, Descriptor>(1, 1, rescaleEngine->clone()),
            cornerNN, coarseLattice, fineLattice, processorLevelCoarse);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXYCorner3D<T, Descriptor>(1, -1, rescaleEngine->clone()),
            cornerNP, coarseLattice, fineLattice, processorLevelCoarse);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXYCorner3D<T, Descriptor>(-1, 1, rescaleEngine->clone()),
            cornerPN, coarseLattice, fineLattice, processorLevelCoarse);

        integrateProcessingFunctional(
            new ScalarCubicInterpolationXYCorner3D<T, Descriptor>(-1, -1, rescaleEngine->clone()),
            cornerPP, coarseLattice, fineLattice, processorLevelCoarse);
    }

    // Add a data processor which imposes a time-cyclic behavior on the fine grid
    //   boundary-dynamics.
    integrateProcessingFunctional(
        new Copy_t1_to_t0_3D<T, Descriptor>(numTimeSteps, executionTime),
        fineGridInterface.multiply(2), fineLattice);

    delete rescaleEngine;
}

template <typename T>
std::vector<MultiScalarField3D<T> *> generateScalarFields(
    MultiGridManagement3D const &management, std::vector<BlockCommunicator3D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics)
{
    std::vector<std::vector<Box3D> > bulks = management.getBulks();
    std::vector<std::vector<plint> > mpiProcesses = management.getMpiProcesses();
    PLB_PRECONDITION(mpiProcesses.size() == 0 || mpiProcesses.size() == bulks.size());
    plint envelopeWidth = 1;
    std::vector<MultiScalarField3D<T> *> scalars(bulks.size());
    for (plint iLevel = 0; iLevel < (plint)bulks.size(); ++iLevel) {
        PLB_ASSERT(mpiProcesses.size() == 0 || mpiProcesses[iLevel].size() == bulks[iLevel].size());
        ExplicitThreadAttribution *threadAttribution = new ExplicitThreadAttribution;
        SparseBlockStructure3D geometry(management.getBoundingBox(iLevel));
        for (pluint iBlock = 0; iBlock < bulks[iLevel].size(); ++iBlock) {
            plint mpiProcess = mpiProcesses.size() > 0 ? mpiProcesses[iLevel][iBlock] : 0;
            plint blockId = geometry.nextIncrementalId();
            geometry.addBlock(
                bulks[iLevel][iBlock], bulks[iLevel][iBlock], blockId);  // TODO:uniqueBulk
            threadAttribution->addBlock(blockId, mpiProcess);
        }
        MultiBlockManagement3D blockManagement(geometry, threadAttribution, envelopeWidth, iLevel);
        scalars[iLevel] = new MultiScalarField3D<T>(
            blockManagement, communicators[iLevel], combinedStatistics[iLevel],
            defaultMultiBlockPolicy3D().getMultiScalarAccess<T>());
    }

    return scalars;
}

template <typename T, int nDim>
std::vector<MultiTensorField3D<T, nDim> *> generateTensorFields(
    MultiGridManagement3D const &management, std::vector<BlockCommunicator3D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics)
{
    std::vector<MultiTensorField3D<T, nDim> *> fields(management.getNumLevels());
    std::vector<std::vector<Box3D> > bulks = management.getBulks();
    std::vector<std::vector<plint> > ids = management.getMpiProcesses();
    PLB_PRECONDITION(ids.size() == 0 || ids.size() == bulks.size());
    plint envelopeWidth = 1;
    fields.resize(bulks.size());
    for (plint iLevel = 0; iLevel < (plint)bulks.size(); ++iLevel) {
        ExplicitThreadAttribution *threadAttribution = new ExplicitThreadAttribution;
        PLB_ASSERT(ids.size() == 0 || ids[iLevel].size() == bulks[iLevel].size());
        SparseBlockStructure3D geometry(management.getBoundingBox(iLevel));
        for (pluint iBlock = 0; iBlock < bulks[iLevel].size(); ++iBlock) {
            plint mpiProcess = ids.size() > 0 ? ids[iLevel][iBlock] : 0;
            plint blockId = geometry.nextIncrementalId();
            geometry.addBlock(bulks[iLevel][iBlock], blockId);
            threadAttribution->addBlock(blockId, mpiProcess);
        }

        MultiBlockManagement3D blockManagement(geometry, threadAttribution, envelopeWidth, iLevel);
        fields[iLevel] = new MultiTensorField3D<T, nDim>(
            blockManagement, communicators[iLevel], combinedStatistics[iLevel],
            defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>());
    }

    return fields;
}

}  // namespace plb

#endif  // MULTI_GRID_GENERATOR_3D_HH
