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
 * Serial implementation of scalar, vector and tensor fields for 2D data analysis.
 * -- header file
 */

#ifndef ATOMIC_CONTAINER_BLOCK_2D_H
#define ATOMIC_CONTAINER_BLOCK_2D_H

#include "atomicBlock/atomicBlock2D.h"
#include "core/globalDefs.h"

namespace plb {

class AtomicContainerDataTransfer2D : public BlockDataTransfer2D {
public:
    AtomicContainerDataTransfer2D() { }
    virtual plint staticCellSize() const
    {
        return 0;
    }
    virtual void send(Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const { }
    virtual void receive(Box2D domain, std::vector<char> const &buffer, modif::ModifT kind) { }
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot2D offset)
    {
        receive(domain, buffer, kind);
    }
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    { }
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind)
    { }
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind,
        Dot2D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, from, kind);
    }
};

class ContainerBlockData {
public:
    ContainerBlockData() : uniqueID(0) { }
    virtual ~ContainerBlockData() { }
    virtual ContainerBlockData *clone() const = 0;
    void setUniqueID(plint uniqueID_)
    {
        uniqueID = uniqueID_;
    }
    plint getUniqueID() const
    {
        return uniqueID;
    }

private:
    plint uniqueID;
};

class AtomicContainerBlock2D : public AtomicBlock2D {
public:
    AtomicContainerBlock2D(plint nx_, plint ny_);
    ~AtomicContainerBlock2D();
    AtomicContainerBlock2D &operator=(AtomicContainerBlock2D const &rhs);
    AtomicContainerBlock2D(AtomicContainerBlock2D const &rhs);
    void swap(AtomicContainerBlock2D &rhs);

public:
    void setData(ContainerBlockData *data_);
    ContainerBlockData *getData();
    ContainerBlockData const *getData() const;
    virtual identifiers::BlockId getBlockId() const;

public:
    virtual BlockDataTransfer2D &getDataTransfer()
    {
        return dataTransfer;
    }
    virtual BlockDataTransfer2D const &getDataTransfer() const
    {
        return dataTransfer;
    }

private:
    ContainerBlockData *data;
    AtomicContainerDataTransfer2D dataTransfer;
};

}  // namespace plb

#endif  // ATOMIC_CONTAINER_BLOCK_2D_H
