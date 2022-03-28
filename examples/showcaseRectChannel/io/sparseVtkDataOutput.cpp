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

#include "io/sparseVtkDataOutput.h"

#include "io/imageWriter.h"
#include "io/parallelIO.h"
#include "io/vtkDataOutput.h"
#include "io/vtkDataOutput.hh"

namespace plb {

SparseVtkImageOutput3D::SparseVtkImageOutput3D(FileName fName_) :
    fName(fName_.defaultPath(global::directories().getVtkOutDir()).defaultExt("vtm")),
    vtmFile(fName.get().c_str())
{
// If we are in a Unix-like or a Windows system we create a directory, called dName, to store the
// atomic block vti files. rdName is the relative path of dName, from where the vtm files are
// stored.
#ifdef PLB_USE_POSIX
    dName = fName.getPath() + "/" + fName.getName();
    rdName = fName.getName();
    bool abortIfExists = false;
    makeDirectory(dName, abortIfExists);
#else
#ifdef PLB_WINDOWS
    dName = fName.getPath() + "\\" + fName.getName();
    rdName = fName.getName();
    bool abortIfExists = false;
    makeDirectory(dName, abortIfExists);
#endif
#endif

    vtmFile << "<?xml version=\"1.0\"?>\n";
#ifdef PLB_BIG_ENDIAN
    vtmFile << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"BigEndian\">\n";
#else
    vtmFile
        << "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
#endif
    vtmFile << "<vtkMultiBlockDataSet>\n";
}

SparseVtkImageOutput3D::~SparseVtkImageOutput3D()
{
    vtmFile << "</vtkMultiBlockDataSet>\n</VTKFile>\n";
}

void SparseVtkImageOutput3D::writeField(
    MultiBlock3D const &multiBlock, plint atomicBlockId, Box3D bulk, std::string fieldName,
    VtkDataWriter3D &vtkOut)
{
    AtomicBlock3D const &atomicBlock = multiBlock.getComponent(atomicBlockId);
    Dot3D location = atomicBlock.getLocation();
    Box3D localBulk = bulk.shift(-location.x, -location.y, -location.z);
    std::vector<std::string> typeInfo = multiBlock.getTypeInfo();
    // Make sure it is a scalar- or tensor-field, but not a lattice.
    PLB_ASSERT(typeInfo.size() == 1);
    std::string scalarType = typeInfo[0];
    std::string blockName = multiBlock.getBlockName();
#ifdef PLB_DEBUG
    bool isTensorFieldName = (blockName.find("TensorField3D_") != std::string::npos);
#endif
    PLB_ASSERT(blockName == "ScalarField3D" || isTensorFieldName || blockName == "NTensorField3D");
    plint nDim = multiBlock.getCellDim();
    if (scalarType == "float") {
        vtkOut.writeDataField<float>(
            atomicBlock.getBlockSerializer(localBulk, IndexOrdering::backward), fieldName, nDim);
    } else if (scalarType == "double") {
        vtkOut.writeDataField<double>(
            atomicBlock.getBlockSerializer(localBulk, IndexOrdering::backward), fieldName, nDim);
    } else if (scalarType == "int") {
        vtkOut.writeDataField<int>(
            atomicBlock.getBlockSerializer(localBulk, IndexOrdering::backward), fieldName, nDim);
    } else {
        pcout << "Error: type \"" << scalarType << "\" not available for VTK output." << std::endl;
        PLB_ASSERT(false);
    }
}

std::string SparseVtkImageOutput3D::writeAtomicBlock(
    Group3D &group, plint partId, plint vtkBlockId, plint atomicBlockId, Box3D bulk, bool pointData,
    bool isLocal, double deltaX, Array<double, 3> offset)
{
    plint numFields = group.getNumBlocks();
    std::string partName(
        fName.getName() + createFileName("-g", vtkBlockId, 3) + createFileName("-", partId, 8));
    FileName fullName(fName);
#ifdef PLB_USE_POSIX
    fullName.setPath(dName);
#else
#ifdef PLB_WINDOWS
    fullName.setPath(dName);
#endif
#endif
    fullName.setName(partName);
    fullName.setExt("vti");
    if (isLocal) {
        bool mainProcOnly = false;
        VtkDataWriter3D vtkOut(fullName, pointData, mainProcOnly);
        if (pointData) {
            Box3D bbox(group.getBoundingBox());
            if (bulk.x1 < bbox.x1) {
                ++bulk.x1;
            }
            if (bulk.y1 < bbox.y1) {
                ++bulk.y1;
            }
            if (bulk.z1 < bbox.z1) {
                ++bulk.z1;
            }
        }
        vtkOut.writeHeader(bulk, offset, deltaX);
        vtkOut.startPiece(bulk);
        for (plint i = 0; i < numFields; ++i) {
            MultiBlock3D const &multiBlock = group.get(i);
            std::string fieldName = group.getName(i);
            writeField(multiBlock, atomicBlockId, bulk, fieldName, vtkOut);
        }
        vtkOut.endPiece();
        vtkOut.writeFooter();
    }
#ifdef PLB_USE_POSIX
    fullName.setPath(rdName);
#else
#ifdef PLB_WINDOWS
    fullName.setPath(rdName);
#else
    fullName.setPath("");
#endif
#endif
    return fullName.get();
}

void SparseVtkImageOutput3D::writeVtkBlock(
    Group3D &group, double deltaX, Array<double, 3> const &offset, plint vtkBlockId, bool pointData)
{
    vtmFile << "<Block index=\"" << vtkBlockId << "\" "
            << "name=\"" << createFileName("Level ", vtkBlockId, 3) << "\">\n";
    MultiBlockManagement3D const &management = group.getMultiBlockManagement();
    ThreadAttribution const &threadAttribution = management.getThreadAttribution();
    SparseBlockStructure3D const &sparseBlock = management.getSparseBlockStructure();
    std::map<plint, Box3D> const &domains = sparseBlock.getBulks();
    std::map<plint, Box3D>::const_iterator it = domains.begin();
    plint partId = 0;
    for (; it != domains.end(); ++it, ++partId) {
        std::string partName = writeAtomicBlock(
            group, partId, vtkBlockId, it->first, it->second, pointData,
            threadAttribution.isLocal(it->first), deltaX, offset);
        vtmFile << "<DataSet index=\"" << partId << "\" "
                << "file=\"" << partName << "\">\n"
                << "</DataSet>\n";
    }
    vtmFile << "</Block>\n";
}

void SparseVtkImageOutput3D::writeVtkBlock(
    Group3D &group, double deltaX, plint vtkBlockId, bool pointData)
{
    Array<double, 3> offset(0.0, 0.0, 0.0);
    writeVtkBlock(group, deltaX, offset, vtkBlockId, pointData);
}

}  // namespace plb
