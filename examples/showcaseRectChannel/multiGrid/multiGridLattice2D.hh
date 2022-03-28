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

#ifndef MULTI_GRID_LATTICE_2D_HH
#define MULTI_GRID_LATTICE_2D_HH

#include "core/util.h"
#include "io/parallelIO.h"
#include "multiGrid/multiGridLattice2D.h"
#include "multiGrid/multiGridParameterManager.h"

namespace plb {

/* ****** MultiGridLattice2D ****** */

template <typename T, template <typename U> class Descriptor>
MultiGridLattice2D<T, Descriptor>::MultiGridLattice2D(
    MultiGridManagement2D management, Dynamics<T, Descriptor> *backgroundDynamics,
    plint behaviorLevel, FineGridInterfaceInstantiator<T, Descriptor> *fineGridInstantiator_,
    CoarseGridInterfaceInstantiator<T, Descriptor> *coarseGridInstantiator_) :
    MultiGrid2D(management, behaviorLevel),
    fineGridInstantiator(fineGridInstantiator_),
    coarseGridInstantiator(coarseGridInstantiator_)
{
    std::vector<Dynamics<T, Descriptor> *> dynamicsVector(this->getNumLevels());
    for (plint iLevel = 0; iLevel < this->getNumLevels(); ++iLevel) {
        Dynamics<T, Descriptor> *cloned = backgroundDynamics->clone();
        cloned->rescale(behaviorLevel - iLevel, behaviorLevel - iLevel);
        dynamicsVector[iLevel] = cloned;
    }
    delete backgroundDynamics;
    // generate the lattices that correspond to the blocks contained inside the management object
    lattices = generateLattices(
        this->getMultiGridManagement(), dynamicsVector,
        defaultMultiGridPolicy2D().getBlockCommunicator<T>(management.getNumLevels()),
        defaultMultiGridPolicy2D().getCombinedStatistics(management.getNumLevels()),
        fineGridInstantiator->getEnvelopeWidth());

    // initialize the internalStatSubscription (multiGrid) once all the lattices have been created
    this->internalStatSubscription().initialize();
}

template <typename T, template <typename U> class Descriptor>
MultiGridLattice2D<T, Descriptor>::MultiGridLattice2D(
    MultiGridManagement2D management, std::vector<BlockCommunicator2D *> communicators_,
    std::vector<CombinedStatistics *> combinedStatistics_,
    Dynamics<T, Descriptor> *backgroundDynamics, plint behaviorLevel,
    FineGridInterfaceInstantiator<T, Descriptor> *fineGridInstantiator_,
    CoarseGridInterfaceInstantiator<T, Descriptor> *coarseGridInstantiator_) :
    MultiGrid2D(management, behaviorLevel),
    fineGridInstantiator(fineGridInstantiator_),
    coarseGridInstantiator(coarseGridInstantiator_)
{
    std::vector<Dynamics<T, Descriptor> *> dynamicsVector(this->getNumLevels());
    for (plint iLevel = 0; iLevel < this->getNumLevels(); ++iLevel) {
        Dynamics<T, Descriptor> *cloned = backgroundDynamics->clone();
        cloned->rescale(behaviorLevel - iLevel, behaviorLevel - iLevel);
        dynamicsVector[iLevel] = cloned;
    }
    delete backgroundDynamics;
    // generate the lattices that correspond to the blocks contained inside the management object
    lattices = generateLattices(
        this->getMultiGridManagement(), dynamicsVector, communicators_, combinedStatistics_,
        fineGridInstantiator->getEnvelopeWidth());
    // initialize the internalStatSubscription (multiGrid) once all the lattices have been created
    this->internalStatSubscription().initialize();
}

/// Get rid properly of all the pointers in the class
template <typename T, template <typename U> class Descriptor>
MultiGridLattice2D<T, Descriptor>::~MultiGridLattice2D()
{
    for (plint i = 0; i < this->getNumLevels(); ++i) {
        delete lattices[i];
        // the dynamics are taken care by each lattice
    }

    delete fineGridInstantiator;
    delete coarseGridInstantiator;
}

template <typename T, template <typename U> class Descriptor>
MultiGridLattice2D<T, Descriptor>::MultiGridLattice2D(
    MultiGridLattice2D<T, Descriptor> const &rhs) :
    MultiGrid2D(rhs),
    fineGridInstantiator(rhs.fineGridInstantiator->clone()),
    coarseGridInstantiator(rhs.coarseGridInstantiator->clone())
{
    std::vector<Dynamics<T, Descriptor> *> dynamicsVector(this->getNumLevels());

    for (plint iLevel = 0; iLevel < this->getNumLevels(); ++iLevel) {
        Dynamics<T, Descriptor> *cloned = rhs.lattices[iLevel]->getBackgroundDynamics().clone();
        dynamicsVector[iLevel] = cloned;
    }
    MultiGridManagement2D const &management = this->getMultiGridManagement();
    // generate the lattices that correspond to the blocks contained inside the management object
    lattices = generateLattices(
        management, dynamicsVector,
        defaultMultiGridPolicy2D().getBlockCommunicator<T>(management.getNumLevels()),
        defaultMultiGridPolicy2D().getCombinedStatistics(management.getNumLevels()),
        fineGridInstantiator->getEnvelopeWidth());

    // initialize the internalStatSubscription (multiGrid) once all the lattices have been created
    this->internalStatSubscription().initialize();
}

template <typename T, template <typename U> class Descriptor>
MultiGridLattice2D<T, Descriptor>::MultiGridLattice2D(
    MultiGridLattice2D<T, Descriptor> const &rhs, Box2D subDomain, bool crop) :
    MultiGrid2D(rhs, subDomain, crop),
    fineGridInstantiator(rhs.fineGridInstantiator),
    coarseGridInstantiator(rhs.coarseGridInstantiator)
{
    std::vector<Dynamics<T, Descriptor> *> backgroundDynamics;
    for (plint iDyn = 0; iDyn < this->getNumLevels(); ++iDyn) {
        backgroundDynamics[iDyn] = new NoDynamics<T, Descriptor>();
    }
    // generate the lattices that correspond to the blocks contained inside the management object
    lattices = generateLattices(
        this->getMultiGridManagement(), backgroundDynamics, this->getCommunicators(),
        this->getCombinedStatistics(), fineGridInstantiator->getEnvelopeWidth());

    // initialize the internalStatSubscription (multiGrid) once all the lattices have been created
    this->internalStatSubscription().initialize();
}

template <typename T, template <typename U> class Descriptor>
MultiGridLattice2D<T, Descriptor>::MultiGridLattice2D(MultiGrid2D const &rhs) : MultiGrid2D(rhs)
{
    std::vector<Dynamics<T, Descriptor> *> backgroundDynamics;
    for (plint iDyn = 0; iDyn < this->getNumLevels(); ++iDyn) {
        backgroundDynamics[iDyn] = new NoDynamics<T, Descriptor>();
    }
    // generate the lattices that correspond to the blocks contained inside the management object
    lattices = generateLattices(
        this->getMultiGridManagement(), backgroundDynamics, this->getCommunicators(),
        this->getCombinedStatistics());

    // initialize the internalStatSubscription (multiGrid) once all the lattices have been created
    this->internalStatSubscription().initialize();
}

template <typename T, template <typename U> class Descriptor>
MultiGridLattice2D<T, Descriptor>::MultiGridLattice2D(
    MultiGrid2D const &rhs, Box2D subDomain, bool crop) :
    MultiGrid2D(rhs, subDomain, crop)
{
    std::vector<Dynamics<T, Descriptor> *> backgroundDynamics;
    for (plint iDyn = 0; iDyn < this->getNumLevels(); ++iDyn) {
        backgroundDynamics[iDyn] = new NoDynamics<T, Descriptor>();
    }
    // generate the lattices that correspond to the blocks contained inside the management object
    lattices = generateLattices(
        this->getMultiGridManagement(), backgroundDynamics, this->getCommunicators(),
        this->getCombinedStatistics());

    // initialize the internalStatSubscription (multiGrid) once all the lattices have been created
    this->internalStatSubscription().initialize();
}

template <typename T, template <typename U> class Descriptor>
MultiGridLattice2D<T, Descriptor> &MultiGridLattice2D<T, Descriptor>::operator=(
    MultiGridLattice2D<T, Descriptor> const &rhs)
{
    MultiGridLattice2D<T, Descriptor> tmp(rhs);
    swap(tmp);
    return *this;
}

/** This method should be called to create the coupling among the internal multi blocks of the
 *  multi grid.
 */
template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::initialize()
{
    this->createInterfaces();

    // eliminate the statistics in the refinement overlap
    this->eliminateStatisticsInOverlap();

    // Execute this data processor once, because this allocates memory on the fine grid dynamcis
    //   for time t0 (after we have values for t1)
    lattices[0]->initialize();
    std::vector<std::vector<Box2D> > fineGridInterfaces =
        this->getMultiGridManagement().getFineInterface();
    plint timeForImmediateExecution = 0;
    plint numTimeSteps = 2;
    // fineGrid interfaces start at level 1
    for (pluint iLevel = 1; iLevel < fineGridInterfaces.size(); ++iLevel) {
        std::vector<Box2D> levelInterfaces = fineGridInterfaces[iLevel];
        for (pluint iInterf = 0; iInterf < levelInterfaces.size(); ++iInterf) {
            applyProcessingFunctional(
                new Copy_t1_to_t0_2D<T, Descriptor>(numTimeSteps, timeForImmediateExecution),
                levelInterfaces[iInterf].multiply(2), *lattices[iLevel]);
        }
        lattices[iLevel]->initialize();
    }

    // toggle all stats on
    for (int iLevel = 0; iLevel < (plint)lattices.size(); ++iLevel) {
        lattices[iLevel]->toggleInternalStatistics(true);
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::createInterfaces()
{
    MultiGridManagement2D management = this->getMultiGridManagement();

    std::vector<std::vector<Box2D> > fineGridInterfaces = management.getFineInterface();
    std::vector<std::vector<Box2D> > coarseGridInterfaces = management.getCoarseInterface();
    std::vector<std::vector<Array<plint, 2> > > coarseOrientations =
        management.getCoarseOrientation();
    std::vector<std::vector<Array<plint, 2> > > fineOrientations = management.getFineOrientation();

    // Create the part of the coupling in which a lattice acts as a coarse lattice.
    for (pluint iCoarse = 0; iCoarse < lattices.size() - 1; ++iCoarse) {
        for (pluint iInterface = 0; iInterface < coarseGridInterfaces[iCoarse].size(); ++iInterface)
        {
            plint direction =
                (coarseOrientations[iCoarse][iInterface])[0];  // w.r.t. the coarse grid
            plint orientation =
                (coarseOrientations[iCoarse][iInterface])[1];  // w.r.t. the coarse grid
            coarseGridInstantiator->instantiateDataProcessors(
                coarseGridInterfaces[iCoarse][iInterface], *lattices[iCoarse],
                *lattices[iCoarse + 1], direction, orientation);
        }
    }
    // Create the part of the coupling in which a lattice acts as a fine lattice.
    for (pluint iFine = 1; iFine < lattices.size(); ++iFine) {
        for (pluint iInterface = 0; iInterface < fineGridInterfaces[iFine].size(); ++iInterface) {
            plint direction = (fineOrientations[iFine][iInterface])[0];    // w.r.t. the fine grid
            plint orientation = (fineOrientations[iFine][iInterface])[1];  // w.r.t. the fine grid
            fineGridInstantiator->instantiateDataProcessors(
                fineGridInterfaces[iFine][iInterface], *lattices[iFine - 1], *lattices[iFine],
                direction, orientation);
        }
    }
}

/* *** MultiGrid2D methods *** */
/// Retrieve the lattices representing each a refinement level
template <typename T, template <typename U> class Descriptor>
MultiBlockLattice2D<T, Descriptor> &MultiGridLattice2D<T, Descriptor>::getComponent(plint iBlock)
{
    PLB_PRECONDITION(iBlock >= 0 && iBlock < (plint)lattices.size());
    return *lattices[iBlock];
}

template <typename T, template <typename U> class Descriptor>
const MultiBlockLattice2D<T, Descriptor> &MultiGridLattice2D<T, Descriptor>::getComponent(
    plint iBlock) const
{
    PLB_PRECONDITION(iBlock >= 0 && iBlock < (plint)lattices.size());
    return *lattices[iBlock];
}

/* **** BlockLatticeBase2D methods **** */

template <typename T, template <typename U> class Descriptor>
Cell<T, Descriptor> &MultiGridLattice2D<T, Descriptor>::get(plint iX, plint iY)
{
    return lattices[this->getBehaviorLevel()]->get(iX, iY);
}

template <typename T, template <typename U> class Descriptor>
Cell<T, Descriptor> const &MultiGridLattice2D<T, Descriptor>::get(plint iX, plint iY) const
{
    return lattices[this->getBehaviorLevel()]->get(iX, iY);
}

template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::specifyStatisticsStatus(Box2D domain, bool status)
{
    for (plint iLevel = 0; iLevel < (plint)lattices.size(); ++iLevel) {
        Box2D rescaledDomain = this->getScaleManager().scaleBox(domain, iLevel);
        lattices[iLevel]->specifyStatisticsStatus(rescaledDomain, status);
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::collide(Box2D domain)
{
    plint levelNumber = (plint)lattices.size();
    PLB_PRECONDITION(levelNumber >= 0 && levelNumber <= (plint)lattices.size());
    for (plint iLevel = 0; iLevel < levelNumber; ++iLevel) {
        // each level will iterate 2 times the iterations of level - 1
        plint fineGridIt = util::roundToInt(util::twoToThePower(iLevel));
        for (plint iterations = 0; iterations < fineGridIt; ++iterations) {
            lattices[iLevel]->collide(domain);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::collide()
{
    plint levelNumber = (plint)lattices.size();
    PLB_PRECONDITION(levelNumber >= 0 && levelNumber <= (plint)lattices.size());
    for (plint iLevel = 0; iLevel < levelNumber; ++iLevel) {
        // each level will iterate 2 times the iterations of level - 1
        plint fineGridIt = util::roundToInt(util::twoToThePower(iLevel));
        for (plint iterations = 0; iterations < fineGridIt; ++iterations) {
            lattices[iLevel]->collide();
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::stream(Box2D domain)
{
    plint levelNumber = (plint)lattices.size();
    PLB_PRECONDITION(levelNumber >= 0 && levelNumber <= (plint)lattices.size());
    for (plint iLevel = 0; iLevel < levelNumber; ++iLevel) {
        // each level will iterate 2 times the iterations of level - 1
        plint fineGridIt = util::roundToInt(util::twoToThePower(iLevel));
        for (plint iterations = 0; iterations < fineGridIt; ++iterations) {
            lattices[iLevel]->stream(domain);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::stream()
{
    plint levelNumber = (plint)lattices.size();
    PLB_PRECONDITION(levelNumber >= 0 && levelNumber <= (plint)lattices.size());
    for (plint iLevel = 0; iLevel < levelNumber; ++iLevel) {
        // each level will iterate 2 times the iterations of level - 1
        plint fineGridIt = util::roundToInt(util::twoToThePower(iLevel));
        for (plint iterations = 0; iterations < fineGridIt; ++iterations) {
            lattices[iLevel]->stream();
        }
        // update the values for parallelism
        lattices[iLevel]->getBlockCommunicator().duplicateOverlaps(
            *lattices[iLevel], modif::staticVariables);
    }

    this->evaluateStatistics();
}

/// One iteration of the multigrid over a given domain (this domain reffers to the reference level)
template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::collideAndStream(Box2D domain)
{
    plint levelNumber = (plint)lattices.size();
    PLB_PRECONDITION(levelNumber >= 0 && levelNumber <= (plint)lattices.size());
    for (plint iLevel = 0; iLevel < levelNumber; ++iLevel) {
        // each level will iterate 2 times the iterations of level - 1
        plint fineGridIt = util::roundToInt(util::twoToThePower(iLevel));
        for (plint iterations = 0; iterations < fineGridIt; ++iterations) {
            lattices[iLevel]->collideAndStream(domain);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::iterateMultiGrid(plint level)
{
    PLB_PRECONDITION(level >= 0 && level < (plint)lattices.size());
    lattices[level]->collideAndStream();
    if ((pluint)level < lattices.size() - 1) {
        iterateMultiGrid(level + 1);
        iterateMultiGrid(level + 1);
        // Overlaps must be duplicated on the coarse lattice, because the
        //   fine->coarse copy acts on bulk nodes only.
        lattices[level]->getBlockCommunicator().duplicateOverlaps(
            *lattices[level], modif::staticVariables);
    }
}

/// One iteration of the entire multigrid
template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::collideAndStream()
{
    iterateMultiGrid(0);
    this->evaluateStatistics();
}

/** This function does not have a sense in the multigrid case. It is
 *  therefore void but implemented to complain with the BlockLatticeBase2D interface
 */
template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::incrementTime()
{ }

template <typename T, template <typename U> class Descriptor>
TimeCounter &MultiGridLattice2D<T, Descriptor>::getTimeCounter()
{
    return lattices[this->getBehaviorLevel()]->getTimeCounter();
}

template <typename T, template <typename U> class Descriptor>
TimeCounter const &MultiGridLattice2D<T, Descriptor>::getTimeCounter() const
{
    return lattices[this->getBehaviorLevel()]->getTimeCounter();
}

template <typename T, template <typename U> class Descriptor>
int MultiGridLattice2D<T, Descriptor>::getBlockId() const
{
    return lattices[this->getReferenceLevel()]->getStaticId();
}

/** Use the information contained in manualGridRefiner to avoid the
 *  computation of statistics in the overlapping regions between the coarse
 *  and fine grids.
 */
template <typename T, template <typename U> class Descriptor>
void MultiGridLattice2D<T, Descriptor>::eliminateStatisticsInOverlap()
{
    // the fine interface is the region inside the fine grid that overlaps with the coarse
    //   grid. It is therefore there that we need to turn off the statistics in the
    //   coarse grid.
    std::vector<std::vector<Box2D> > fineOverlaps =
        this->getMultiGridManagement().getFineInterface();
    for (plint iLevel = 0; iLevel < (plint)lattices.size() - 1; ++iLevel) {
        for (plint iOv = 0; iOv < (plint)fineOverlaps[iLevel].size(); ++iOv) {
            Box2D currentOverlap = (fineOverlaps[iLevel])[iOv];
            lattices[iLevel]->specifyStatisticsStatus(currentOverlap, false);
        }
    }
    // it is also necessary to avoid computation in the coarse interface, as the refinement
    //   overlap is two coarse sites
    std::vector<std::vector<Box2D> > coarseOverlaps =
        this->getMultiGridManagement().getCoarseInterface();
    for (plint iLevel = 0; iLevel < (plint)lattices.size() - 1; ++iLevel) {
        for (plint iOv = 0; iOv < (plint)coarseOverlaps[iLevel].size(); ++iOv) {
            Box2D currentOverlap = (coarseOverlaps[iLevel])[iOv];
            lattices[iLevel]->specifyStatisticsStatus(currentOverlap, false);
        }
    }
}

/** Interpolate and decimate the multi blocks that form the MultiGridLattice2D
 *   in order to have several multi blocks of the same level. Then join them by
 *   performing a union of multi blocks.
 */
template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >
    MultiGridLattice2D<T, Descriptor>::convertToLevel(plint level) const
{
    // create the resulting multiBlock
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > result =
        generateMultiBlockLattice<T, Descriptor>(
            lattices[level]->getBoundingBox(), lattices[level]->getBackgroundDynamics().clone(),
            lattices[level]->getMultiBlockManagement().getEnvelopeWidth());

    // create the first lattice to start looping
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > refined =
        generateMultiBlockLattice<T, Descriptor>(
            lattices[0]->getBoundingBox(), lattices[0]->getBackgroundDynamics().clone(),
            lattices[0]->getMultiBlockManagement().getEnvelopeWidth());

    copyNonLocal<T, Descriptor>(
        *lattices[0], *refined, refined->getBoundingBox(), modif::staticVariables);

    // loop to interpolate the blocks until level-1
    for (plint iLevel = 0; iLevel < (plint)level; ++iLevel) {
        plint envelopeWidth = lattices[iLevel + 1]->getMultiBlockManagement().getEnvelopeWidth();
        // interpolate lattice at iLevel
        std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > tmp =
            refine(*refined, -1, -1, refined->getBackgroundDynamics().clone());

        refined = generateMultiBlockLattice<T, Descriptor>(
            lattices[iLevel + 1]->getBoundingBox(),
            lattices[iLevel + 1]->getBackgroundDynamics().clone(), envelopeWidth);

        // copy from reshaped iLevel to result
        copyNonLocal<T, Descriptor>(
            *tmp, *refined, refined->getBoundingBox(), modif::staticVariables);
        // copy from iLevel+1 to result
        copyNonLocal<T, Descriptor>(
            *lattices[iLevel + 1], *refined, refined->getBoundingBox(), modif::staticVariables);
    }

    // create the last lattice to start looping in the other orientation
    plint lastLevel = getNumLevels() - 1;
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > coarsened =
        std::unique_ptr<MultiBlockLattice2D<T, Descriptor> >(
            new MultiBlockLattice2D<T, Descriptor>(*lattices[lastLevel]));

    defineDynamics<T, Descriptor>(
        *coarsened, coarsened->getBoundingBox(),
        lattices[lastLevel]->getBackgroundDynamics().clone());

    copyNonLocal<T, Descriptor>(
        *lattices[lastLevel], *coarsened, coarsened->getBoundingBox(), modif::staticVariables);

    // loop to decimate the blocks until level+1
    for (plint iLevel = getNumLevels() - 1; iLevel >= (plint)level + 1; --iLevel) {
        // interpolate lattice at iLevel
        std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > tmp =
            coarsen(*coarsened, 1, 1, lattices[iLevel]->getBackgroundDynamics().clone());

        coarsened = generateJoinMultiBlockLattice<T, Descriptor>(*tmp, *lattices[iLevel - 1]);
        // make all the dynamics from refinement disappear
        defineDynamics<T, Descriptor>(
            *coarsened, coarsened->getBoundingBox(),
            lattices[iLevel - 1]->getBackgroundDynamics().clone());

        // copy from reshaped iLevel to result
        copyNonLocal<T, Descriptor>(
            *tmp, *coarsened, coarsened->getBoundingBox(), modif::staticVariables);
        // copy from iLevel-1 to result
        copyNonLocal<T, Descriptor>(
            *lattices[iLevel - 1], *coarsened, coarsened->getBoundingBox(), modif::staticVariables);
    }

    // final copies
    copyNonLocal<T, Descriptor>(
        *refined, *result, result->getBoundingBox(), modif::staticVariables);
    copyNonLocal<T, Descriptor>(
        *lattices[level], *result, result->getBoundingBox(), modif::staticVariables);
    copyNonLocal<T, Descriptor>(
        *coarsened, *result, result->getBoundingBox(), modif::staticVariables);

    return result;
}

/////////// Free Functions //////////////////////////////

template <typename T, template <typename U> class Descriptor>
double getStoredAverageDensity(MultiGridLattice2D<T, Descriptor> const &multiGrid)
{
    return Descriptor<T>::fullRho(
        multiGrid.getInternalStatistics().getAverage(LatticeStatistics::avRhoBar));
}

template <typename T, template <typename U> class Descriptor>
double getStoredAverageEnergy(MultiGridLattice2D<T, Descriptor> const &multiGrid)
{
    return 0.5 * multiGrid.getInternalStatistics().getAverage(LatticeStatistics::avUSqr);
}

template <typename T, template <typename U> class Descriptor>
double getStoredMaxVelocity(MultiGridLattice2D<T, Descriptor> const &multiGrid)
{
    return std::sqrt(multiGrid.getInternalStatistics().getMax(LatticeStatistics::maxUSqr));
}

}  // namespace plb

#endif  // MULTI_GRID_LATTICE_2D_HH
