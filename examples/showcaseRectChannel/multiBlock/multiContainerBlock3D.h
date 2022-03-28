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

#ifndef MULTI_CONTAINER_BLOCK_3D_H
#define MULTI_CONTAINER_BLOCK_3D_H

#include "atomicBlock/atomicContainerBlock3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"

namespace plb {

class MultiContainerBlock3D : public MultiBlock3D {
public:
    typedef std::map<plint, AtomicContainerBlock3D *> BlockMap;

public:
    MultiContainerBlock3D(
        MultiBlockManagement3D const &multiBlockManagement_,
        CombinedStatistics *combinedStatistics_);
    MultiContainerBlock3D(plint nx_, plint ny_, plint nz_);
    MultiContainerBlock3D(MultiBlock3D const &rhs);
    MultiContainerBlock3D(MultiBlock3D const &rhs, Box3D subDomain, bool crop);
    ~MultiContainerBlock3D();
    MultiContainerBlock3D &operator=(MultiContainerBlock3D const &rhs);
    MultiContainerBlock3D(MultiContainerBlock3D const &rhs);
    MultiContainerBlock3D *clone() const;
    MultiContainerBlock3D *clone(MultiBlockManagement3D const &multiBlockManagement) const;
    void swap(MultiContainerBlock3D &rhs);

public:
    virtual AtomicContainerBlock3D &getComponent(plint iBlock);
    virtual AtomicContainerBlock3D const &getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);
    std::string getBlockName() const;
    std::vector<std::string> getTypeInfo() const;

private:
    void allocateBlocks();
    void allocateBlocks(MultiContainerBlock3D const &rhs);
    void deAllocateBlocks();

private:
    BlockMap blocks;
};

MultiContainerBlock3D *createContainerBlock(MultiBlock3D &templ, ContainerBlockData *data);
MultiContainerBlock3D *createContainerBlock(
    MultiBlock3D &templ, ContainerBlockData *data, plint envelopeWidth);

}  // namespace plb

#endif  // MULTI_CONTAINER_BLOCK_3D_H
