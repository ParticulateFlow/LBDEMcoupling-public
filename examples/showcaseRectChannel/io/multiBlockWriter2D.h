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
 * I/O routines for 2D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_WRITER_2D_H
#define MULTI_BLOCK_WRITER_2D_H

#include "atomicBlock/atomicContainerBlock2D.h"
#include "core/globalDefs.h"
#include "core/serializer.h"
#include "io/plbFiles.h"
#include "multiBlock/multiBlock2D.h"

namespace plb {

namespace parallelIO {

void save(MultiBlock2D &multiBlock, FileName fName, bool dynamicContent = true);

void saveFull(
    MultiBlock2D &multiBlock, FileName fName, IndexOrdering::OrderingT = IndexOrdering::forward);

/** Convert the multi-block into raw serialized data, but preserve the
 *  parallel distribution. This is a preparation for parallel IO. The
 *  IDs of the original blocks are not preserved; they are renumbered
 *  contiguously, in increasing order of the original IDs.
 *  @var offset: Position, in bytes, at which the data corresponding
 *               to a given atomic-block must be written into the file.
 *               The size of the vector equals the number of atomic-blocks.
 *  @var myBlockIds: New, contiguously numbered IDs of the atomic-blocks
 *                   which are local to the current MPI thread.
 *  @var data: The serialized data of the local atomic-blocks.
 **/
void dumpData(
    MultiBlock2D &multiBlock, bool dynamicContent, std::vector<plint> &offset,
    std::vector<plint> &myBlockIds, std::vector<std::vector<char> > &data);

void writeXmlSpec(
    MultiBlock2D &multiBlock, FileName fName, std::vector<plint> const &offset,
    bool dynamicContent);

}  // namespace parallelIO

}  // namespace plb

#endif  // MULTI_BLOCK_WRITER_2D_H
