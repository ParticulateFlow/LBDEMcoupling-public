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

#ifndef MULTI_GRID_LATTICE_3D_H
#define MULTI_GRID_LATTICE_3D_H

#include <memory>

#include "core/array.h"
#include "core/block3D.h"
#include "core/blockStatistics.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "multiBlock/blockCommunicator3D.h"
#include "multiBlock/combinedStatistics.h"
#include "multiBlock/multiBlock3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiGrid/gridRefinement.h"
#include "multiGrid/gridRefinementDynamics.h"
#include "multiGrid/multiGrid3D.h"
#include "multiGrid/multiGridParameterManager.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class MultiGridLattice3D : public BlockLatticeBase3D<T, Descriptor>, public MultiGrid3D {
public:
    MultiGridLattice3D(
        MultiGridManagement3D management, std::vector<BlockCommunicator3D *> communicators_,
        std::vector<CombinedStatistics *> combinedStatistics_,
        Dynamics<T, Descriptor> *backgroundDynamics, plint behaviorLevel = 0);

    MultiGridLattice3D(
        MultiGridManagement3D management, Dynamics<T, Descriptor> *backgroundDynamics,
        plint behaviorLevel = 0);

    /// Copy constructor for the whole multi grid
    MultiGridLattice3D(MultiGridLattice3D<T, Descriptor> const &rhs);
    /// Copy constructor for a subdomain of the multi grid
    MultiGridLattice3D(
        MultiGridLattice3D<T, Descriptor> const &rhs, Box3D subDomain, bool crop = true);

    MultiGridLattice3D(MultiGrid3D const &rhs);
    MultiGridLattice3D(MultiGrid3D const &rhs, Box3D subDomain, bool crop = true);

    MultiGridLattice3D<T, Descriptor> &operator=(MultiGridLattice3D<T, Descriptor> const &rhs);

    ~MultiGridLattice3D();

    /// Create the couplings between lattices
    void initialize();

    void createInterfaces();

    /// Create a single multiBlock that represents the multiGrid.
    std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > convertToLevel(plint level) const;

    /* *** MultiGrid3D methods *** */
    int getBlockId() const;

    /// Retrieve the lattices representing each a refinement level
    MultiBlockLattice3D<T, Descriptor> &getComponent(plint iBlock);
    const MultiBlockLattice3D<T, Descriptor> &getComponent(plint iBlock) const;

    /* **** BlockLatticeBase3D methods **** */

    virtual Cell<T, Descriptor> &get(plint iX, plint iY, plint iZ);
    virtual Cell<T, Descriptor> const &get(plint iX, plint iY, plint iZ) const;
    virtual void specifyStatisticsStatus(Box3D domain, bool status);
    virtual void collide(Box3D domain);
    virtual void collide();
    virtual void stream(Box3D domain);
    virtual void stream();
    virtual void collideAndStream(Box3D domain);
    virtual void collideAndStream();
    virtual void incrementTime();
    TimeCounter &getTimeCounter();
    TimeCounter const &getTimeCounter() const;

private:
    void iterateMultiGrid(plint level);
    void eliminateStatisticsInOverlap();

private:
    std::vector<MultiBlockLattice3D<T, Descriptor> *> lattices;
};

template <typename T, template <typename U> class Descriptor>
double getStoredAverageDensity(MultiGridLattice3D<T, Descriptor> const &multiGrid);

template <typename T, template <typename U> class Descriptor>
double getStoredAverageEnergy(MultiGridLattice3D<T, Descriptor> const &multiGrid);

template <typename T, template <typename U> class Descriptor>
double getStoredMaxVelocity(MultiGridLattice3D<T, Descriptor> const &multiGrid);

}  // namespace plb

#endif  // MULTI_GRID_LATTICE_3D_H
