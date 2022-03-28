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

#ifndef MULTI_CONTAINER_BLOCK_2D_H
#define MULTI_CONTAINER_BLOCK_2D_H

#include "atomicBlock/atomicContainerBlock2D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlock2D.h"

namespace plb {

class MultiContainerBlock2D : public MultiBlock2D {
public:
    typedef std::map<plint, AtomicContainerBlock2D *> BlockMap;

public:
    MultiContainerBlock2D(
        MultiBlockManagement2D const &multiBlockManagement_,
        CombinedStatistics *combinedStatistics_);
    MultiContainerBlock2D(plint nx_, plint ny);
    MultiContainerBlock2D(MultiBlock2D const &rhs);
    MultiContainerBlock2D(MultiBlock2D const &rhs, Box2D subDomain, bool crop);
    ~MultiContainerBlock2D();
    MultiContainerBlock2D &operator=(MultiContainerBlock2D const &rhs);
    MultiContainerBlock2D(MultiContainerBlock2D const &rhs);
    MultiContainerBlock2D *clone() const;
    MultiContainerBlock2D *clone(MultiBlockManagement2D const &multiBlockManagement) const;
    void swap(MultiContainerBlock2D &rhs);
    std::string getBlockName() const;
    std::vector<std::string> getTypeInfo() const;

public:
    virtual AtomicContainerBlock2D &getComponent(plint iBlock);
    virtual AtomicContainerBlock2D const &getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);

private:
    void allocateBlocks();
    void allocateBlocks(MultiContainerBlock2D const &rhs);
    void deAllocateBlocks();

private:
    BlockMap blocks;
};

MultiContainerBlock2D *createContainerBlock(MultiBlock2D &templ, ContainerBlockData *data);

}  // namespace plb

#endif  // MULTI_CONTAINER_BLOCK_2D_H
