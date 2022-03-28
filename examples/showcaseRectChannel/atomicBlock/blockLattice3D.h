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
 * The dynamics of a 3D block lattice -- header file.
 */
#ifndef BLOCK_LATTICE_3D_H
#define BLOCK_LATTICE_3D_H

#include <map>
#include <vector>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/blockIdentifiers.h"
#include "core/blockLatticeBase3D.h"
#include "core/cell.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"

/// All Palabos code is contained in this namespace.
namespace plb {

template <typename T, template <typename U> class Descriptor>
struct Dynamics;
template <typename T, template <typename U> class Descriptor>
class BlockLattice3D;

template <typename T, template <typename U> class Descriptor>
class BlockLatticeDataTransfer3D : public BlockDataTransfer3D {
public:
    BlockLatticeDataTransfer3D();
    virtual void setBlock(AtomicBlock3D &block);
    virtual void setConstBlock(AtomicBlock3D const &block);
    virtual BlockLatticeDataTransfer3D<T, Descriptor> *clone() const;
    virtual plint staticCellSize() const;
    /// Send data from the lattice into a byte-stream.
    virtual void send(Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the lattice.
    virtual void receive(Box3D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset)
    {
        receive(domain, buffer, kind);
    }
    /// Receive data from a byte-stream into the block, and re-map IDs for dynamics if exist.
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds);
    /// Attribute data between two lattices.
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind);
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind, Dot3D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
    }

private:
    void send_static(Box3D domain, std::vector<char> &buffer) const;
    void send_dynamic(Box3D domain, std::vector<char> &buffer) const;
    void send_all(Box3D domain, std::vector<char> &buffer) const;

    void receive_static(Box3D domain, std::vector<char> const &buffer);
    void receive_dynamic(Box3D domain, std::vector<char> const &buffer);
    void receive_all(Box3D domain, std::vector<char> const &buffer);
    void receive_regenerate(
        Box3D domain, std::vector<char> const &buffer,
        std::map<int, int> const &idIndirect = (std::map<int, int>()));

    void attribute_static(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        BlockLattice3D<T, Descriptor> const &from);
    void attribute_dynamic(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        BlockLattice3D<T, Descriptor> const &from);
    void attribute_all(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        BlockLattice3D<T, Descriptor> const &from);
    void attribute_regenerate(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
        BlockLattice3D<T, Descriptor> const &from);

private:
    BlockLattice3D<T, Descriptor> *lattice;
    BlockLattice3D<T, Descriptor> const *constLattice;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class ExternalRhoJcollideAndStream3D;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class ExternalCollideAndStream3D;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class WaveAbsorptionExternalRhoJcollideAndStream3D;
};

/// A regular lattice for highly efficient 3D LB dynamics.
/** A block lattice contains a regular array of Cell objects and
 * some useful methods to execute the LB dynamics on the lattice.
 *
 * This class is not intended to be derived from.
 */
template <typename T, template <typename U> class Descriptor>
class BlockLattice3D : public BlockLatticeBase3D<T, Descriptor>, public AtomicBlock3D {
public:
    /// Construction of an nx_ by ny_ by nz_ lattice
    BlockLattice3D(plint nx_, plint ny_, plint nz_, Dynamics<T, Descriptor> *backgroundDynamics_);
    /// Destruction of the lattice
    ~BlockLattice3D();
    /// Copy construction
    BlockLattice3D(BlockLattice3D<T, Descriptor> const &rhs);
    /// Copy assignment
    BlockLattice3D &operator=(BlockLattice3D<T, Descriptor> const &rhs);
    /// Swap the content of two BlockLattices
    void swap(BlockLattice3D &rhs);

public:
    /// Read/write access to lattice cells
    virtual Cell<T, Descriptor> &get(plint iX, plint iY, plint iZ)
    {
        PLB_PRECONDITION(iX < this->getNx());
        PLB_PRECONDITION(iY < this->getNy());
        PLB_PRECONDITION(iZ < this->getNz());
        return grid[iX][iY][iZ];
    }
    /// Read only access to lattice cells
    virtual Cell<T, Descriptor> const &get(plint iX, plint iY, plint iZ) const
    {
        PLB_PRECONDITION(iX < this->getNx());
        PLB_PRECONDITION(iY < this->getNy());
        PLB_PRECONDITION(iZ < this->getNz());
        return grid[iX][iY][iZ];
    }
    /// Specify wheter statistics measurements are done on a rect. domain
    virtual void specifyStatisticsStatus(Box3D domain, bool status);
    /// Apply collision step to a 3D sub-box
    virtual void collide(Box3D domain);
    /// Apply collision step to the whole domain
    virtual void collide();
    /// Apply streaming step to a 3D sub-box
    virtual void stream(Box3D domain);
    /// Apply streaming step to the whole domain
    virtual void stream();
    /// Apply first collision, then streaming step to a 3D sub-box
    virtual void collideAndStream(Box3D domain);
    /// Apply first collision, then streaming step to the whole domain
    virtual void collideAndStream();
    /// Increment time counter
    /** Warning: don't call this method manually. Instead, call incrementTime()
     *  on the multi-block lattice. Otherwise, the internal time of the multi-block
     *  and the atomic-blocks get out of sync.
     **/
    virtual void incrementTime();

public:
    /// Attribute dynamics to a cell.
    void attributeDynamics(plint iX, plint iY, plint iZ, Dynamics<T, Descriptor> *dynamics);
    /// Get a reference to the background dynamics
    Dynamics<T, Descriptor> &getBackgroundDynamics();
    /// Get a const reference to the background dynamics
    Dynamics<T, Descriptor> const &getBackgroundDynamics() const;
    /// Assign an individual clone of the new dynamics to every cell.
    void resetDynamics(Dynamics<T, Descriptor> const &dynamics);
    /// Apply streaming step to bulk (non-boundary) cells
    void bulkStream(Box3D domain);
    /// Apply streaming step to boundary cells
    void boundaryStream(Box3D bound, Box3D domain);
    /// Apply collision and streaming step to bulk (non-boundary) cells
    void bulkCollideAndStream(Box3D domain);

private:
    /// Generic implementation of bulkCollideAndStream(domain).
    void linearBulkCollideAndStream(Box3D domain);
    /// Cache-efficient implementation of bulkCollideAndStream(domain)for
    ///   nearest-neighbor lattices.
    void blockwiseBulkCollideAndStream(Box3D domain);

private:
    /// Helper method for memory allocation
    void allocateAndInitialize();
    /// Helper method for memory de-allocation
    void releaseMemory();
    void implementPeriodicity();
    plint allocatedMemory() const;

private:
    void periodicDomain(Box3D domain);

private:
    Dynamics<T, Descriptor> *backgroundDynamics;
    Cell<T, Descriptor> *rawData;
    Cell<T, Descriptor> ***grid;

public:
    static CachePolicy3D &cachePolicy();
    template <typename T_, template <typename U_> class Descriptor_>
    friend class ExternalRhoJcollideAndStream3D;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class ExternalCollideAndStream3D;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class PackedExternalRhoJcollideAndStream3D;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class WaveAbsorptionExternalRhoJcollideAndStream3D;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class OnLinkExternalRhoJcollideAndStream3D;
    template <typename T_, template <typename U_> class Descriptor_>
    friend class GrayScaleCollideAndStream3D;
};

template <typename T, template <typename U> class Descriptor>
double getStoredAverageDensity(BlockLattice3D<T, Descriptor> const &blockLattice);

template <typename T, template <typename U> class Descriptor>
double getStoredAverageEnergy(BlockLattice3D<T, Descriptor> const &blockLattice);

template <typename T, template <typename U> class Descriptor>
double getStoredAverageVelocity(BlockLattice3D<T, Descriptor> const &blockLattice);

}  // namespace plb

#endif  // BLOCK_LATTICE_3D_H
