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
 * Geometry specifications for 3D multiblocks -- header file.
 */

#include "multiBlock/nonLocalTransfer3D.h"

#include "core/globalDefs.h"

namespace plb {

std::vector<Overlap3D> copyAllDataTransfer(
    SparseBlockStructure3D const &block1, SparseBlockStructure3D const &block2)
{
    std::vector<Overlap3D> dataTransfer;  // The return value.
    std::vector<plint> ids;               // Temporary value.
    std::vector<Box3D> intersections;     // Temporary value.
    std::map<plint, Box3D>::const_iterator it1 = block1.getBulks().begin();
    for (; it1 != block1.getBulks().end(); ++it1) {
        ids.clear();
        intersections.clear();
        block2.intersect(it1->second, ids, intersections);
        for (pluint iInters = 0; iInters < intersections.size(); ++iInters) {
            dataTransfer.push_back(Overlap3D(it1->first, ids[iInters], intersections[iInters]));
        }
    }
    return dataTransfer;
}

std::vector<Overlap3D> copyDomainDataTransfer(
    SparseBlockStructure3D const &block1, Box3D block1Domain, SparseBlockStructure3D const &block2,
    Box3D block2Domain)
{
    PLB_PRECONDITION(block1Domain.getNx() == block2Domain.getNx());
    PLB_PRECONDITION(block1Domain.getNy() == block2Domain.getNy());
    PLB_PRECONDITION(block1Domain.getNz() == block2Domain.getNz());
    plint shiftX = block1Domain.x0 - block2Domain.x0;
    plint shiftY = block1Domain.y0 - block2Domain.y0;
    plint shiftZ = block1Domain.z0 - block2Domain.z0;
    std::vector<plint> block1Ids, block2Ids;
    std::vector<Box3D> block1Components, block2Inters;
    std::vector<Overlap3D> dataTransfer;  // The return value.
    block1.intersect(block1Domain, block1Ids, block1Components);
    PLB_ASSERT(block1Ids.size() == block1Components.size());
    for (pluint iComp1 = 0; iComp1 < block1Ids.size(); ++iComp1) {
        block2Ids.clear();
        block2Inters.clear();
        block2.intersect(
            block1Components[iComp1].shift(-shiftX, -shiftY, -shiftZ), block2Ids, block2Inters);
        for (pluint iInters = 0; iInters < block2Inters.size(); ++iInters) {
            dataTransfer.push_back(Overlap3D(
                block1Ids[iComp1], block2Ids[iInters],
                block2Inters[iInters].shift(shiftX, shiftY, shiftZ), shiftX, shiftY, shiftZ));
        }
    }

    return dataTransfer;
}

void copy_generic(
    MultiBlock3D const &from, Box3D const &fromDomain, MultiBlock3D &to, Box3D const &toDomain,
    modif::ModifT typeOfModif)
{
    Box3D fromDomain_(fromDomain);
    Box3D toDomain_(toDomain);
    adjustEqualSize(fromDomain_, toDomain_);
    std::vector<Overlap3D> dataTransfer = copyDomainDataTransfer(
        from.getMultiBlockManagement().getSparseBlockStructure(), fromDomain_,
        to.getMultiBlockManagement().getSparseBlockStructure(), toDomain_);
    to.getBlockCommunicator().communicate(dataTransfer, from, to, typeOfModif);
    to.getBlockCommunicator().duplicateOverlaps(to, typeOfModif);
}

// TODO: This one is suspicious, domain is ununsed.
void copyNonLocal_generic(
    MultiBlock3D const &from, MultiBlock3D &to, Box3D const &domain, modif::ModifT typeOfModif)
{
    std::vector<Overlap3D> dataTransfer = copyAllDataTransfer(
        from.getMultiBlockManagement().getSparseBlockStructure(),
        to.getMultiBlockManagement().getSparseBlockStructure());
    to.getBlockCommunicator().communicate(dataTransfer, from, to, typeOfModif);
    to.getBlockCommunicator().duplicateOverlaps(to, typeOfModif);
}

}  // namespace plb
