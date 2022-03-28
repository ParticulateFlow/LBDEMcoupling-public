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

#ifndef MULTI_GRID_LATTICE_2D_H
#define MULTI_GRID_LATTICE_2D_H

#include <memory>

#include "core/array.h"
#include "core/block2D.h"
#include "core/blockStatistics.h"
#include "core/geometry2D.h"
#include "core/globalDefs.h"
#include "multiBlock/blockCommunicator2D.h"
#include "multiBlock/combinedStatistics.h"
#include "multiBlock/multiBlock2D.h"
#include "multiBlock/multiBlockManagement2D.h"
#include "multiGrid/coarseGridProcessors2D.h"
#include "multiGrid/fineGridProcessors2D.h"
#include "multiGrid/gridRefinement.h"
#include "multiGrid/gridRefinementDynamics.h"
#include "multiGrid/multiGrid2D.h"
#include "multiGrid/multiGridGenerator2D.h"
#include "multiGrid/multiGridParameterManager.h"

namespace plb {

/// Main class when dealing with grid refinement
template <typename T, template <typename U> class Descriptor>
class MultiGridLattice2D : public BlockLatticeBase2D<T, Descriptor>, public MultiGrid2D {
public:
    MultiGridLattice2D(
        MultiGridManagement2D management, std::vector<BlockCommunicator2D *> communicators_,
        std::vector<CombinedStatistics *> combinedStatistics_,
        Dynamics<T, Descriptor> *backgroundDynamics, plint behaviorLevel,
        FineGridInterfaceInstantiator<T, Descriptor> *fineGridInstantiator_,
        CoarseGridInterfaceInstantiator<T, Descriptor> *coarseGridInstantiator_);

    MultiGridLattice2D(
        MultiGridManagement2D management, Dynamics<T, Descriptor> *backgroundDynamics,
        plint behaviorLevel, FineGridInterfaceInstantiator<T, Descriptor> *fineGridInstantiator_,
        CoarseGridInterfaceInstantiator<T, Descriptor> *coarseGridInstantiator_);

    /// Copy constructor for the whole multi grid
    MultiGridLattice2D(MultiGridLattice2D<T, Descriptor> const &rhs);
    /// Copy constructor for a subdomain of the multi grid
    MultiGridLattice2D(
        MultiGridLattice2D<T, Descriptor> const &rhs, Box2D subDomain, bool crop = true);

    MultiGridLattice2D(MultiGrid2D const &rhs);
    MultiGridLattice2D(MultiGrid2D const &rhs, Box2D subDomain, bool crop = true);

    MultiGridLattice2D<T, Descriptor> &operator=(MultiGridLattice2D<T, Descriptor> const &rhs);
    ~MultiGridLattice2D();

    /// Create the couplings between lattices
    void initialize();

    /// Create a single multiBlock that represents the multiGrid.
    std::unique_ptr<MultiBlockLattice2D<T, Descriptor> > convertToLevel(plint level) const;

    /* *** MultiGrid2D methods *** */
    int getBlockId() const;

    /// Retrieve the lattices representing each a refinement level
    MultiBlockLattice2D<T, Descriptor> &getComponent(plint iBlock);
    const MultiBlockLattice2D<T, Descriptor> &getComponent(plint iBlock) const;

    /* **** BlockLatticeBase2D methods **** */

    virtual Cell<T, Descriptor> &get(plint iX, plint iY);
    virtual Cell<T, Descriptor> const &get(plint iX, plint iY) const;
    virtual void specifyStatisticsStatus(Box2D domain, bool status);
    virtual void collide(Box2D domain);
    virtual void collide();
    virtual void stream(Box2D domain);
    virtual void stream();
    virtual void collideAndStream(Box2D domain);
    virtual void collideAndStream();
    virtual void incrementTime();
    TimeCounter &getTimeCounter();
    TimeCounter const &getTimeCounter() const;

private:
    void createInterfaces();
    void iterateMultiGrid(plint level);
    void eliminateStatisticsInOverlap();

private:
    std::vector<MultiBlockLattice2D<T, Descriptor> *> lattices;
    // the objects to create the interfaces
    FineGridInterfaceInstantiator<T, Descriptor> *fineGridInstantiator;
    CoarseGridInterfaceInstantiator<T, Descriptor> *coarseGridInstantiator;
};

template <typename T, template <typename U> class Descriptor>
double getStoredAverageDensity(MultiGridLattice2D<T, Descriptor> const &multiGrid);

template <typename T, template <typename U> class Descriptor>
double getStoredAverageEnergy(MultiGridLattice2D<T, Descriptor> const &multiGrid);

template <typename T, template <typename U> class Descriptor>
double getStoredMaxVelocity(MultiGridLattice2D<T, Descriptor> const &multiGrid);

}  // namespace plb

#endif
