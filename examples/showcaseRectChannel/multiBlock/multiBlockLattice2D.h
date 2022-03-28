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
 * A 2D multiblock lattice -- header file.
 */
#ifndef MULTI_BLOCK_LATTICE_2D_H
#define MULTI_BLOCK_LATTICE_2D_H

#include <vector>

#include "atomicBlock/blockLattice2D.h"
#include "core/blockLatticeBase2D.h"
#include "core/blockStatistics.h"
#include "core/cell.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlock2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class BlockLattice2D;

template <typename T, template <typename U> class Descriptor>
struct MultiCellAccess2D {
    virtual ~MultiCellAccess2D() { }
    virtual Cell<T, Descriptor> &getDistributedCell(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, BlockLattice2D<T, Descriptor> *> &lattices) = 0;
    virtual Cell<T, Descriptor> const &getDistributedCell(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, BlockLattice2D<T, Descriptor> *> const &lattices) const = 0;
    virtual void broadCastCell(
        Cell<T, Descriptor> &cell, plint fromBlock,
        MultiBlockManagement2D const &multiBlockManagement) const = 0;
    virtual MultiCellAccess2D<T, Descriptor> *clone() const = 0;
};

/// A complex LatticeBase, itself decomposed into smaller components.
/** This extensible class can be used for example for cache-optimized
 * lattices, irregular domains (no memory allocation in areas exterior to
 * the domain) and parallel lattices. The actual behavior of the lattice
 * is parametrizable by a multiBlockHandler instance, which is given to
 * the constructor.
 *
 * The MultiBlockLattice does not itself possess LatticeProcessors. The Lattice-
 * Processors are delegated to the respective LatticeBases.
 */
template <typename T, template <typename U> class Descriptor>
class MultiBlockLattice2D : public BlockLatticeBase2D<T, Descriptor>, public MultiBlock2D {
public:
    typedef std::map<plint, BlockLattice2D<T, Descriptor> *> BlockMap;

public:
    MultiBlockLattice2D(
        MultiBlockManagement2D const &multiBlockManagement_,
        BlockCommunicator2D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiCellAccess2D<T, Descriptor> *multiCellAccess_,
        Dynamics<T, Descriptor> *backgroundDynamics_);
    MultiBlockLattice2D(plint nx, plint ny, Dynamics<T, Descriptor> *backgroundDynamics_);
    ~MultiBlockLattice2D();
    MultiBlockLattice2D(MultiBlockLattice2D<T, Descriptor> const &rhs);
    MultiBlockLattice2D(MultiBlock2D const &rhs);
    virtual MultiBlockLattice2D<T, Descriptor> *clone() const;
    virtual MultiBlockLattice2D<T, Descriptor> *clone(
        MultiBlockManagement2D const &newManagement) const;
    /// Extract sub-domain from rhs and construct a multi-block-lattice with the same
    ///  data distribution and policy-classes; but the data itself and the data-processors
    ///  are not copied. MultiCellAccess takes default value.
    MultiBlockLattice2D(MultiBlock2D const &rhs, Box2D subDomain, bool crop = true);
    /// Attention: data-processors of rhs, which were pointing at rhs, will continue pointing
    /// to rhs, and not to *this.
    void swap(MultiBlockLattice2D &rhs);
    /// Attention: data-processors of rhs, which were pointing at rhs, will continue pointing
    /// to rhs, and not to *this.
    MultiBlockLattice2D<T, Descriptor> &operator=(MultiBlockLattice2D<T, Descriptor> const &rhs);
    // Assign an individual clone of the new dynamics to every cell.
    void resetDynamics(Dynamics<T, Descriptor> const &dynamics);

    Dynamics<T, Descriptor> const &getBackgroundDynamics() const;
    virtual Cell<T, Descriptor> &get(plint iX, plint iY);
    virtual Cell<T, Descriptor> const &get(plint iX, plint iY) const;
    virtual void specifyStatisticsStatus(Box2D domain, bool status);
    virtual void collide(Box2D domain);
    virtual void collide();
    virtual void stream(Box2D domain);
    virtual void stream();
    virtual void externalStream();
    virtual void collideAndStream(Box2D domain);
    virtual void collideAndStream();
    void externalCollideAndStream();
    virtual void incrementTime();
    virtual void resetTime(pluint value);
    virtual BlockLattice2D<T, Descriptor> &getComponent(plint blockId);
    virtual BlockLattice2D<T, Descriptor> const &getComponent(plint blockId) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);

public:
    BlockMap &getBlockLattices();
    BlockMap const &getBlockLattices() const;
    virtual void getDynamicsDict(Box2D domain, std::map<std::string, int> &dict);
    virtual std::string getBlockName() const;
    virtual std::vector<std::string> getTypeInfo() const;
    static std::string blockName();
    static std::string basicType();
    static std::string descriptorType();

private:
    void collideAndStreamImplementation();
    void streamImplementation();
    void allocateAndInitialize();
    void eliminateStatisticsInEnvelope();
    Box2D extendPeriodic(Box2D const &box, plint envelopeWidth) const;

private:
    Dynamics<T, Descriptor> *backgroundDynamics;
    MultiCellAccess2D<T, Descriptor> *multiCellAccess;
    BlockMap blockLattices;

public:
    static const int staticId;
};

template <typename T, template <typename U> class Descriptor>
MultiBlockLattice2D<T, Descriptor> &findMultiBlockLattice2D(id_t id);

template <typename T, template <typename U> class Descriptor>
double getStoredAverageDensity(MultiBlockLattice2D<T, Descriptor> const &blockLattice);

template <typename T, template <typename U> class Descriptor>
double getStoredAverageEnergy(MultiBlockLattice2D<T, Descriptor> const &blockLattice);

template <typename T, template <typename U> class Descriptor>
double getStoredMaxVelocity(MultiBlockLattice2D<T, Descriptor> const &blockLattice);

}  // namespace plb

#endif  // MULTI_BLOCK_LATTICE_2D_H
