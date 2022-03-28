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

#ifndef MULTI_PARTICLE_FIELD_3D_H
#define MULTI_PARTICLE_FIELD_3D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"
#include "particles/particleField3D.h"

namespace plb {

template <class ParticleFieldT>
class MultiParticleField3D : public MultiBlock3D {
public:
    typedef std::map<plint, ParticleFieldT *> BlockMap;

public:
    MultiParticleField3D(
        MultiBlockManagement3D const &multiBlockManagement_,
        CombinedStatistics *combinedStatistics_);
    MultiParticleField3D(plint nx_, plint ny_, plint nz_);
    MultiParticleField3D(MultiBlock3D const &rhs);
    MultiParticleField3D(MultiBlock3D const &rhs, Box3D subDomain, bool crop);
    ~MultiParticleField3D();
    virtual MultiParticleField3D<ParticleFieldT> *clone() const;
    virtual MultiParticleField3D<ParticleFieldT> *clone(
        MultiBlockManagement3D const &newManagement) const;
    MultiParticleField3D &operator=(MultiParticleField3D<ParticleFieldT> const &rhs);
    MultiParticleField3D(MultiParticleField3D<ParticleFieldT> const &rhs);
    void swap(MultiParticleField3D<ParticleFieldT> &rhs);

public:
    virtual ParticleFieldT &getComponent(plint iBlock);
    virtual ParticleFieldT const &getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);
    std::string getBlockName() const;
    std::vector<std::string> getTypeInfo() const;
    static std::string blockName();
    static std::string basicType();
    static std::string descriptorType();

private:
    void allocateBlocks();
    void deAllocateBlocks();

private:
    BlockMap blocks;

public:
    static const int staticId;
};

template <class ParticleFieldT>
MultiParticleField3D<ParticleFieldT> &findMultiParticleField3D(id_t id);

}  // namespace plb

#endif  // MULTI_PARTICLE_FIELD_3D_H
