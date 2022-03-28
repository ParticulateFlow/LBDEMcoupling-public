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
 * I/O routines for 3D multiblock -- implementation.
 */

#include "io/multiBlockWriter3D.h"

#include <algorithm>
#include <memory>
#include <numeric>

#include "core/globalDefs.h"
#include "core/multiBlockIdentifiers3D.h"
#include "core/plbProfiler.h"
#include "core/plbTypenames.h"
#include "core/processorIdentifiers3D.h"
#include "core/util.h"
#include "io/mpiParallelIO.h"
#include "io/plbFiles.h"
#include "libraryInterfaces/TINYXML_xmlIO.h"
#include "libraryInterfaces/TINYXML_xmlIO.hh"
#include "multiBlock/multiBlockOperations3D.h"
#include "multiBlock/nonLocalTransfer3D.h"
#include "parallelism/mpiManager.h"
#ifdef HDF5
#include "io/hdfWrapper.h"
#endif
#include <ctime>

namespace plb {

namespace parallelIO {

/***** 1. Multi-Block Writer **************************************************/

XMLwriter writeXmlSpec(
    MultiBlock3D &multiBlock, FileName fName, std::vector<plint> const &offset, bool dynamicContent)
{
    MultiBlockManagement3D const &management = multiBlock.getMultiBlockManagement();
    std::map<plint, Box3D> const &bulks = management.getSparseBlockStructure().getBulks();
    PLB_ASSERT(offset.empty() || bulks.size() == offset.size());
    std::vector<std::string> typeInfo = multiBlock.getTypeInfo();
    std::string blockName = multiBlock.getBlockName();
    PLB_ASSERT(!typeInfo.empty());

    XMLwriter xml;
    XMLwriter &xmlMultiBlock = xml["Block3D"];
    xmlMultiBlock["General"]["Family"].setString(blockName);
    xmlMultiBlock["General"]["Datatype"].setString(typeInfo[0]);
    if (typeInfo.size() > 1) {
        xmlMultiBlock["General"]["Descriptor"].setString(typeInfo[1]);
    }
    xmlMultiBlock["General"]["cellDim"].set(multiBlock.getCellDim());
    xmlMultiBlock["General"]["dynamicContent"].set(dynamicContent);
    xmlMultiBlock["General"]["globalId"].set(multiBlock.getId());

    Array<plint, 6> boundingBox = multiBlock.getBoundingBox().to_plbArray();
    xmlMultiBlock["Structure"]["BoundingBox"].set<plint, 6>(boundingBox);
    xmlMultiBlock["Structure"]["EnvelopeWidth"].set(management.getEnvelopeWidth());
    xmlMultiBlock["Structure"]["GridLevel"].set(management.getRefinementLevel());
    xmlMultiBlock["Structure"]["NumComponents"].set(bulks.size());

    xmlMultiBlock["Data"]["File"].setString(FileName(fName).setExt("dat"));

    XMLwriter &xmlBulks = xmlMultiBlock["Data"]["Component"];
    std::map<plint, Box3D>::const_iterator it = bulks.begin();
    plint iComp = 0;
    for (; it != bulks.end(); ++it) {
        Box3D bulk = it->second;
        xmlBulks[iComp].set<plint, 6>(bulk.to_plbArray());
        ++iComp;
    }
    if (!offset.empty()) {
        xmlMultiBlock["Data"]["Offsets"].set(offset);
    }

    // The following prints a unique list of dynamics-id pairs for all dynamics
    //   classes used in the multi-block. This is necessary, because dynamics
    //   classes may be ordered differently from one compilation to the other,
    //   or from one compiler to the other.
    //
    // Given that the dynamics classes are unique, they can be indexed by their
    //   name (which is not the case of the data processors below).
    std::map<std::string, int> dynamicsDict;
    multiBlock.getDynamicsDict(multiBlock.getBoundingBox(), dynamicsDict);
    if (!dynamicsDict.empty()) {
        XMLwriter &xmlDynamicsDict = xmlMultiBlock["Data"]["DynamicsDict"];
        for (std::map<std::string, int>::const_iterator it = dynamicsDict.begin();
             it != dynamicsDict.end(); ++it)
        {
            xmlDynamicsDict[it->first].set(it->second);
        }
    }

    // This is the only section in which actual content is stored outside the
    //   binary blob: the serialization of the data processors. This
    //   serialization was chosen to be in ASCII, because it takes little space
    //   and can be somewhat complicated.
    //
    // It is important that the processors are indexed by a continuous index
    //   "iProcessor". They cannot be indexed by the class name ("Name") or static
    //   id ("id") because several instances of the same class may occur.
    XMLwriter &xmlProcessors = xmlMultiBlock["Data"]["Processor"];
    std::vector<MultiBlock3D::ProcessorStorage3D> const &processors =
        multiBlock.getStoredProcessors();
    for (plint iProcessor = 0; iProcessor < (plint)processors.size(); ++iProcessor) {
        int id = processors[iProcessor].getGenerator().getStaticId();
        if (id >= 0) {
            Box3D domain;
            std::string data;
            processors[iProcessor].getGenerator().serialize(domain, data);
            xmlProcessors[iProcessor]["Name"].set(meta::processorRegistration3D().getName(id));
            xmlProcessors[iProcessor]["Domain"].set<plint, 6>(domain.to_plbArray());
            xmlProcessors[iProcessor]["Data"].setString(data);
            xmlProcessors[iProcessor]["Level"].set(processors[iProcessor].getLevel());
            xmlProcessors[iProcessor]["Blocks"].set(processors[iProcessor].getMultiBlockIds());
        }
    }
    return xml;
}

void writeOneBlockXmlSpec(
    MultiBlock3D const &multiBlock, FileName fName, plint dataSize,
    IndexOrdering::OrderingT ordering)
{
    fName.setExt("plb");
    MultiBlockManagement3D const &management = multiBlock.getMultiBlockManagement();
    std::vector<std::string> typeInfo = multiBlock.getTypeInfo();
    std::string blockName = multiBlock.getBlockName();
    PLB_ASSERT(!typeInfo.empty());

    XMLwriter xml;
    XMLwriter &xmlMultiBlock = xml["Block3D"];
    xmlMultiBlock["General"]["Family"].setString(blockName);
    xmlMultiBlock["General"]["Datatype"].setString(typeInfo[0]);
    if (typeInfo.size() > 1) {
        xmlMultiBlock["General"]["Descriptor"].setString(typeInfo[1]);
    }
    xmlMultiBlock["General"]["cellDim"].set(multiBlock.getCellDim());
    bool dynamicContent = false;
    xmlMultiBlock["General"]["dynamicContent"].set(dynamicContent);

    Array<plint, 6> boundingBox = multiBlock.getBoundingBox().to_plbArray();
    xmlMultiBlock["Structure"]["BoundingBox"].set<plint, 6>(boundingBox);
    xmlMultiBlock["Structure"]["EnvelopeWidth"].set(management.getEnvelopeWidth());
    xmlMultiBlock["Structure"]["NumComponents"].set(1);

    xmlMultiBlock["Data"]["File"].setString(FileName(fName).setExt("dat"));
    if (ordering == IndexOrdering::forward) {
        xmlMultiBlock["Data"]["IndexOrdering"].setString("zIsFastest");
    } else {
        xmlMultiBlock["Data"]["IndexOrdering"].setString("xIsFastest");
    }
    XMLwriter &xmlBulks = xmlMultiBlock["Data"]["Component"];
    xmlBulks.set<plint, 6>(multiBlock.getBoundingBox().to_plbArray());

    xmlMultiBlock["Data"]["Offsets"].set(dataSize);

    xml.print(FileName(fName).defaultPath(global::directories().getOutputDir()));
}

void transposeToBackward(plint sizeOfCell, Box3D const &domain, std::vector<char> &data)
{
    plint nx = domain.getNx();
    plint ny = domain.getNy();
    plint nz = domain.getNz();
    std::vector<char> transp(data.size());

    for (plint iX = 0; iX < nx; ++iX) {
        for (plint iY = 0; iY < ny; ++iY) {
            for (plint iZ = 0; iZ < nz; ++iZ) {
                plint iForward = sizeOfCell * (iZ + nz * (iY + ny * iX));
                plint iBackward = sizeOfCell * (iX + nx * (iY + ny * iZ));
                for (plint iByte = 0; iByte < sizeOfCell; ++iByte) {
                    transp[iBackward + iByte] = data[iForward + iByte];
                }
            }
        }
    }
    transp.swap(data);
}

#ifdef HDF5
void saveHDF(MultiBlock3D &multiBlock, FileName fName, bool dynamicContent)
{
    global::profiler().start("io");
    std::vector<plint> offset;
    std::vector<plint> myBlockIds;
    std::vector<std::vector<char>> data;

    dumpData(multiBlock, dynamicContent, offset, myBlockIds, data);

    XMLwriter xml = writeXmlSpec(multiBlock, fName, offset, dynamicContent);
    std::string s_fname = fName.get();
    std::string s_xml = xml.sprint();

    writeParallelHDF5(
        s_fname.c_str(), "binary_blob", myBlockIds, offset, data, s_xml, global::mpi().getRank(),
        global::mpi().getGlobalCommunicator());

    global::profiler().stop("io");
}
#endif

void save(MultiBlock3D &multiBlock, FileName fName, bool dynamicContent)
{
    global::profiler().start("io");
    std::vector<plint> offset;
    std::vector<plint> myBlockIds;
    std::vector<std::vector<char>> data;

    dumpData(multiBlock, dynamicContent, offset, myBlockIds, data);

    writeRawData(fName, myBlockIds, offset, data);

    XMLwriter xml = writeXmlSpec(multiBlock, fName, offset, dynamicContent);
    fName.setExt("plb");
    xml.print(FileName(fName).defaultPath(global::directories().getOutputDir()));
    global::profiler().stop("io");
}

void saveFull(
    MultiBlock3D &multiBlock, FileName fName, IndexOrdering::OrderingT ordering, bool appendMode)
{
    global::profiler().start("io");
    SparseBlockStructure3D blockStructure(multiBlock.getBoundingBox());
    Box3D bbox = multiBlock.getBoundingBox();
    if (ordering == IndexOrdering::forward) {
        plint nBlocks = std::min(bbox.getNx(), (plint)global::mpi().getSize());
        std::vector<std::pair<plint, plint>> ranges;
        util::linearRepartition(bbox.x0, bbox.x1, nBlocks, ranges);
        for (pluint iRange = 0; iRange < ranges.size(); ++iRange) {
            blockStructure.addBlock(
                Box3D(
                    ranges[iRange].first, ranges[iRange].second, bbox.y0, bbox.y1, bbox.z0,
                    bbox.z1),
                iRange);
        }
    } else if (ordering == IndexOrdering::backward) {
        plint nBlocks = std::min(bbox.getNz(), (plint)global::mpi().getSize());
        std::vector<std::pair<plint, plint>> ranges;
        util::linearRepartition(bbox.z0, bbox.z1, nBlocks, ranges);
        for (pluint iRange = 0; iRange < ranges.size(); ++iRange) {
            blockStructure.addBlock(
                Box3D(
                    bbox.x0, bbox.x1, bbox.y0, bbox.y1, ranges[iRange].first,
                    ranges[iRange].second),
                iRange);
        }
    } else {
        // Sparse ordering not defined.
        PLB_ASSERT(false);
    }
    plint envelopeWidth = 1;
    MultiBlockManagement3D adjacentMultiBlockManagement(
        blockStructure, new OneToOneThreadAttribution, envelopeWidth);
    MultiBlock3D *multiAdjacentBlock = multiBlock.clone(adjacentMultiBlockManagement);

    std::vector<plint> offset;
    std::vector<plint> myBlockIds;
    std::vector<std::vector<char>> data;

    bool dynamicContent = false;
    dumpData(*multiAdjacentBlock, dynamicContent, offset, myBlockIds, data);
    if (ordering == IndexOrdering::backward && myBlockIds.size() == 1) {
        PLB_ASSERT(data.size() == 1);
        Box3D domain;
        blockStructure.getBulk(myBlockIds[0], domain);
        plint sizeOfCell = multiAdjacentBlock->sizeOfCell();
        PLB_ASSERT(domain.nCells() * sizeOfCell == (plint)data[0].size());
        transposeToBackward(sizeOfCell, domain, data[0]);
    }

    plint totalSize = offset[offset.size() - 1];
    if (!appendMode) {
        writeOneBlockXmlSpec(*multiAdjacentBlock, fName, totalSize, ordering);
    }

    writeRawData(fName, myBlockIds, offset, data, appendMode);
    delete multiAdjacentBlock;
    global::profiler().stop("io");
}

void dumpData(
    MultiBlock3D &multiBlock, bool dynamicContent, std::vector<plint> &offset,
    std::vector<plint> &myBlockIds, std::vector<std::vector<char>> &data)
{
    MultiBlockManagement3D const &management = multiBlock.getMultiBlockManagement();
    std::map<plint, Box3D> const &bulks = management.getSparseBlockStructure().getBulks();

    plint numBlocks = (plint)bulks.size();
    std::map<plint, plint> toContiguousId;
    std::map<plint, Box3D>::const_iterator it = bulks.begin();
    plint pos = 0;
    for (; it != bulks.end(); ++it) {
        toContiguousId[it->first] = pos;
        ++pos;
    }

    std::vector<plint> const &myBlocks = management.getLocalInfo().getBlocks();
    myBlockIds.resize(myBlocks.size());
    data.resize(myBlocks.size());
    std::vector<plint> blockSize(numBlocks);
    std::fill(blockSize.begin(), blockSize.end(), 0);

    for (pluint iBlock = 0; iBlock < myBlocks.size(); ++iBlock) {
        plint blockId = myBlocks[iBlock];
        SmartBulk3D bulk(management, blockId);
        Box3D localBulk(bulk.toLocal(bulk.getBulk()));
        // localBulk is the real domain without the communication buffer
        AtomicBlock3D const &block = multiBlock.getComponent(blockId);
        modif::ModifT typeOfVariables =
            dynamicContent ? modif::dataStructure : modif::staticVariables;
        block.getDataTransfer().send(localBulk, data[iBlock], typeOfVariables);
        plint contiguousId = toContiguousId[blockId];
        myBlockIds[iBlock] = contiguousId;
        blockSize[contiguousId] = (plint)data[iBlock].size();
    }
#ifdef PLB_MPI_PARALLEL
    global::mpi().allReduceVect(blockSize, MPI_SUM);
#endif
    offset.resize(numBlocks);
    std::partial_sum(blockSize.begin(), blockSize.end(), offset.begin());
}

}  // namespace parallelIO

}  // namespace plb
