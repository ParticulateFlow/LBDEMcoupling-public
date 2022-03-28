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
 * Serial implementation of scalar, vector and tensor fields for 3D data analysis.
 * -- header file
 */

#ifndef ATOMIC_CONTAINER_BLOCK_3D_H
#define ATOMIC_CONTAINER_BLOCK_3D_H

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/atomicContainerBlock2D.h"
#include "core/globalDefs.h"

namespace plb {

class AtomicContainerDataTransfer3D : public BlockDataTransfer3D {
public:
    AtomicContainerDataTransfer3D() { }
    virtual void setBlock(AtomicBlock3D &block_) { }
    virtual void setConstBlock(AtomicBlock3D const &block_) { }
    virtual AtomicContainerDataTransfer3D *clone() const
    {
        return new AtomicContainerDataTransfer3D(*this);
    }
    virtual plint staticCellSize() const
    {
        return 0;
    }
    virtual void send(Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const { }
    virtual void receive(Box3D domain, std::vector<char> const &buffer, modif::ModifT kind) { }
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    { }
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset)
    {
        receive(domain, buffer, kind);
    }
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind)
    { }
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind, Dot3D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
    }
};

class AtomicContainerBlock3D : public AtomicBlock3D {
public:
    AtomicContainerBlock3D(plint nx_, plint ny_, plint nz_);
    ~AtomicContainerBlock3D();
    AtomicContainerBlock3D &operator=(AtomicContainerBlock3D const &rhs);
    AtomicContainerBlock3D(AtomicContainerBlock3D const &rhs);
    void swap(AtomicContainerBlock3D &rhs);

public:
    void setData(ContainerBlockData *data_);
    ContainerBlockData *getData();
    ContainerBlockData const *getData() const;
    virtual identifiers::BlockId getBlockId() const;

private:
    ContainerBlockData *data;
};

}  // namespace plb

#endif  // ATOMIC_CONTAINER_BLOCK_3D_H
