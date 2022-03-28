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

#ifndef VTK_DATA_OUTPUT_HH
#define VTK_DATA_OUTPUT_HH

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <typeinfo>

#include "core/globalDefs.h"
#include "core/util.h"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "dataProcessors/ntensorAnalysisWrapper2D.h"
#include "dataProcessors/ntensorAnalysisWrapper2D.hh"
#include "dataProcessors/ntensorAnalysisWrapper3D.h"
#include "io/plbFiles.h"
#include "io/serializerIO.h"
#include "io/vtkDataOutput.h"
#include "parallelism/mpiManager.h"

namespace plb {

////////// struct VtkTypeNames ///////////////////////////////////

template <typename T>
class VtkTypeNames {
public:
    static std::string getName();

private:
    static std::string getBaseName();
};

template <typename T>
std::string VtkTypeNames<T>::getName()
{
    std::stringstream sstream;
    sstream << getBaseName();
    sstream << 8 * sizeof(T);

    std::string tName;
    sstream >> tName;
    return tName;
}

////////// class VtkDataWriter3D ////////////////////////////////////////

template <typename T>
void VtkDataWriter3D::writeDataField(
    DataSerializer const *serializer, std::string const &name, plint nDim)
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        (*ostr) << "<DataArray type=\"" << VtkTypeNames<T>::getName() << "\" Name=\"" << name
                << "\" format=\"binary\" encoding=\"base64";
        if (nDim > 1) {
            (*ostr) << "\" NumberOfComponents=\"" << nDim;
        }
        (*ostr) << "\">\n";
    }

    // Undocumented requirement of the vtk xml file format:
    // in front of every binary blob, base64 or raw-binary, appended or not,
    // there is an UInt32 length indicator, giving the size of the binary blob in bytes;
    // when using base64 encoding, that length header must be encoded separately;
    // there must be no newline between the encoded length indicator and the encoded data block.

    bool enforceUint =
        true;  // VTK uses "unsigned" to indicate the size of data, even on a 64-bit machine.
    serializerToBase64Stream(serializer, ostr, enforceUint, mainProcOnly);

    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        (*ostr) << "\n</DataArray>\n";
    }
}

////////// class VtkDataWriter3D ////////////////////////////////////////

template <typename T>
void VtkAsciiDataWriter3D::writeDataField(Block3D const &block, std::string const &name, plint nDim)
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        (*ostr) << "<DataArray type=\"" << VtkTypeNames<T>::getName() << "\" Name=\"" << name
                << "\" format=\"ascii";
        if (nDim > 1) {
            (*ostr) << "\" NumberOfComponents=\"" << nDim;
        }
        (*ostr) << "\">\n";
    }

    plint numDigits = 0;  // Number of digits is handled by iostream manipulators,
                          // as opposed to being imposed here.
    serializerToAsciiStream<T>(
        block.getBlockSerializer(block.getBoundingBox(), IndexOrdering::backward), ostr, numDigits);

    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        (*ostr) << "\n</DataArray>\n";
    }
}

////////// class ParallelVtkDataWriter3D ////////////////////////////////////////

template <typename T>
void ParallelVtkDataWriter3D::declareDataField(
    std::string const &name, plint nDim, plint declarationPos, plint dataOffset)
{
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(
            fileName.c_str(), std::ios_base::in | std::ios_base::out | std::ios_base::binary);
        ostr->seekp(declarationPos);
        (*ostr) << "\n<DataArray type=\"" << VtkTypeNames<T>::getName() << "\" Name=\"" << name
                << "\" format=\"appended";
        if (nDim > 1) {
            (*ostr) << "\" NumberOfComponents=\"" << nDim;
        }
        (*ostr) << "\" offset=\"" << dataOffset << "\"/>\n";
        delete ostr;
    }
}

////////// class VtkImageOutput2D ////////////////////////////////////

template <typename T>
VtkImageOutput2D<T>::VtkImageOutput2D(std::string fName, double deltaX_) :
    fullName(global::directories().getVtkOutDir() + fName + ".vti"),
    vtkOut(fullName),
    deltaX(deltaX_),
    offset(T(), T()),
    headerWritten(false)
{ }

template <typename T>
VtkImageOutput2D<T>::VtkImageOutput2D(std::string fName, double deltaX_, Array<double, 2> offset_) :
    fullName(global::directories().getVtkOutDir() + fName + ".vti"),
    vtkOut(fullName),
    deltaX(deltaX_),
    offset(offset_),
    headerWritten(false)
{ }

template <typename T>
VtkImageOutput2D<T>::~VtkImageOutput2D()
{
    writeFooter();
}

template <typename T>
void VtkImageOutput2D<T>::writeHeader(plint nx_, plint ny_)
{
    if (headerWritten) {
        PLB_PRECONDITION(nx == nx_);
        PLB_PRECONDITION(ny == ny_);
    } else {
        nx = nx_;
        ny = ny_;
        vtkOut.writeHeader(
            Box3D(0, nx - 1, 0, ny - 1, 0, 0), Array<T, 3>(offset[0], offset[1], T()), deltaX);
        vtkOut.startPiece(Box3D(0, nx - 1, 0, ny - 1, 0, 0));
        headerWritten = true;
    }
}

template <typename T>
void VtkImageOutput2D<T>::writeHeader(Box2D boundingBox_)
{
    Box3D boundingBox3D(boundingBox_.x0, boundingBox_.x1, boundingBox_.y0, boundingBox_.y1, 0, 0);
    if (headerWritten) {
        PLB_PRECONDITION(boundingBox == boundingBox3D);
    } else {
        boundingBox = boundingBox3D;
        vtkOut.writeHeader(boundingBox3D, Array<T, 3>(offset[0], offset[1], T()), deltaX);
        vtkOut.startPiece(boundingBox3D);
        headerWritten = true;
    }
}

template <typename T>
void VtkImageOutput2D<T>::writeFooter()
{
    if (headerWritten) {
        vtkOut.endPiece();
        vtkOut.writeFooter();
        headerWritten = false;
    }
}

template <typename T>
template <typename TConv>
void VtkImageOutput2D<T>::writeData(
    plint nx, plint ny, plint nDim, DataSerializer const *serializer, std::string const &name)
{
    writeHeader(nx, ny);
    vtkOut.writeDataField<TConv>(serializer, name, nDim);
}

template <typename T>
template <typename TConv>
void VtkImageOutput2D<T>::writeData(
    Box2D boundingBox, plint nDim, DataSerializer const *serializer, std::string const &name)
{
    writeHeader(boundingBox);
    vtkOut.writeDataField<TConv>(serializer, name, nDim);
}

template <typename T>
template <typename TConv>
void VtkImageOutput2D<T>::writeData(
    ScalarField2D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    std::unique_ptr<ScalarField2D<TConv> > transformedField = copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    writeData<TConv>(
        scalarField.getNx(), scalarField.getNy(), 1,
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        scalarFieldName);
}

template <typename T>
template <typename TConv>
void VtkImageOutput2D<T>::writeData(
    MultiScalarField2D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    std::unique_ptr<MultiScalarField2D<TConv> > transformedField =
        copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    writeData<TConv>(
        scalarField.getBoundingBox(), 1,
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        scalarFieldName);
}

template <typename T>
template <plint n, typename TConv>
void VtkImageOutput2D<T>::writeData(
    TensorField2D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    std::unique_ptr<TensorField2D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    writeData<TConv>(
        tensorField.getNx(), tensorField.getNy(), n,
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        tensorFieldName);
}

template <typename T>
template <plint n, typename TConv>
void VtkImageOutput2D<T>::writeData(
    MultiTensorField2D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    std::unique_ptr<MultiTensorField2D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    writeData<TConv>(
        tensorField.getBoundingBox(), n,
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        tensorFieldName);
}

template <typename T>
template <typename TConv>
void VtkImageOutput2D<T>::writeData(
    MultiNTensorField2D<T> &nTensorField, std::string nTensorFieldName)
{
    MultiNTensorField2D<TConv> *transformedField =
        copyConvert<T, TConv>(nTensorField, nTensorField.getBoundingBox());
    writeData<TConv>(
        transformedField->getNx(), transformedField->getNy(), transformedField->getNdim(),
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        nTensorFieldName);
    delete transformedField;
}

////////// class VtkImageOutput3D ////////////////////////////////////

template <typename T>
VtkImageOutput3D<T>::VtkImageOutput3D(std::string fName, double deltaX_) :
    fullName(FileName(fName + ".vti").defaultPath(global::directories().getVtkOutDir())),
    vtkOut(fullName),
    deltaX(deltaX_),
    offset(T(), T(), T()),
    headerWritten(false)
{ }

template <typename T>
VtkImageOutput3D<T>::VtkImageOutput3D(std::string fName, double deltaX_, Array<double, 3> offset_) :
    fullName(FileName(fName + ".vti").defaultPath(global::directories().getVtkOutDir())),
    vtkOut(fullName),
    deltaX(deltaX_),
    offset(offset_),
    headerWritten(false)
{ }

template <typename T>
VtkImageOutput3D<T>::~VtkImageOutput3D()
{
    writeFooter();
}

template <typename T>
void VtkImageOutput3D<T>::writeHeader(plint nx_, plint ny_, plint nz_)
{
    writeHeader(Box3D(0, nx_ - 1, 0, ny_ - 1, 0, nz_ - 1));
}

template <typename T>
void VtkImageOutput3D<T>::writeHeader(Box3D boundingBox_)
{
    if (headerWritten) {
        PLB_PRECONDITION(boundingBox == boundingBox_);
    } else {
        boundingBox = boundingBox_;
        vtkOut.writeHeader(boundingBox, offset, deltaX);
        vtkOut.startPiece(boundingBox);
        headerWritten = true;
    }
}

template <typename T>
void VtkImageOutput3D<T>::writeFooter()
{
    if (headerWritten) {
        vtkOut.endPiece();
        vtkOut.writeFooter();
        headerWritten = false;
    }
}

template <typename T>
template <typename TConv>
void VtkImageOutput3D<T>::writeData(
    plint nx, plint ny, plint nz, plint nDim, DataSerializer const *serializer,
    std::string const &name)
{
    writeHeader(nx, ny, nz);
    vtkOut.writeDataField<TConv>(serializer, name, nDim);
}

template <typename T>
template <typename TConv>
void VtkImageOutput3D<T>::writeData(
    Box3D boundingBox, plint nDim, DataSerializer const *serializer, std::string const &name)
{
    writeHeader(boundingBox);
    vtkOut.writeDataField<TConv>(serializer, name, nDim);
}

template <typename T>
template <typename TConv>
void VtkImageOutput3D<T>::writeData(
    ScalarField3D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    std::unique_ptr<ScalarField3D<TConv> > transformedField = copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    writeData<TConv>(
        scalarField.getNx(), scalarField.getNy(), scalarField.getNz(), 1,
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        scalarFieldName);
}

template <typename T>
template <typename TConv>
void VtkImageOutput3D<T>::writeData(
    MultiScalarField3D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    std::unique_ptr<MultiScalarField3D<TConv> > transformedField =
        copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    writeData<TConv>(
        scalarField.getBoundingBox(), 1,
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        scalarFieldName);
}

template <typename T>
template <plint n, typename TConv>
void VtkImageOutput3D<T>::writeData(
    TensorField3D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    std::unique_ptr<TensorField3D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    writeData<TConv>(
        tensorField.getNx(), tensorField.getNy(), tensorField.getNz(), n,
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        tensorFieldName);
}

template <typename T>
template <plint n, typename TConv>
void VtkImageOutput3D<T>::writeData(
    MultiTensorField3D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    std::unique_ptr<MultiTensorField3D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    writeData<TConv>(
        tensorField.getBoundingBox(), n,
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        tensorFieldName);
}

template <typename T>
template <typename TConv>
void VtkImageOutput3D<T>::writeData(
    MultiNTensorField3D<T> &nTensorField, std::string nTensorFieldName)
{
    MultiNTensorField3D<TConv> *transformedField =
        copyConvert<T, TConv>(nTensorField, nTensorField.getBoundingBox());
    writeData<TConv>(
        transformedField->getBoundingBox(), transformedField->getNdim(),
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        nTensorFieldName);
    delete transformedField;
}

////////// class VtkAsciiImageOutput3D ////////////////////////////////////

template <typename T>
VtkAsciiImageOutput3D<T>::VtkAsciiImageOutput3D(
    std::string fName, std::string const &scalars_, double deltaX_) :
    fullName(global::directories().getVtkOutDir() + fName + ".vti"),
    vtkOut(fullName),
    scalars(scalars_),
    deltaX(deltaX_),
    offset(T(), T(), T()),
    headerWritten(false)
{ }

template <typename T>
VtkAsciiImageOutput3D<T>::VtkAsciiImageOutput3D(
    std::string fName, std::string const &scalars_, double deltaX_, Array<double, 3> offset_) :
    fullName(global::directories().getVtkOutDir() + fName + ".vti"),
    vtkOut(fullName),
    scalars(scalars_),
    deltaX(deltaX_),
    offset(offset_),
    headerWritten(false)
{ }

template <typename T>
VtkAsciiImageOutput3D<T>::~VtkAsciiImageOutput3D()
{
    writeFooter();
}

template <typename T>
void VtkAsciiImageOutput3D<T>::writeHeader(plint nx_, plint ny_, plint nz_)
{
    writeHeader(Box3D(0, nx_ - 1, 0, ny_ - 1, 0, nz_ - 1));
}

template <typename T>
void VtkAsciiImageOutput3D<T>::writeHeader(Box3D boundingBox_)
{
    if (headerWritten) {
        PLB_PRECONDITION(boundingBox == boundingBox_);
    } else {
        boundingBox = boundingBox_;
        vtkOut.writeHeader(boundingBox, offset, deltaX);
        vtkOut.startPiece(boundingBox, scalars);
        headerWritten = true;
    }
}

template <typename T>
void VtkAsciiImageOutput3D<T>::writeFooter()
{
    if (headerWritten) {
        vtkOut.endPiece();
        vtkOut.writeFooter();
        headerWritten = false;
    }
}

template <typename T>
template <typename TConv>
void VtkAsciiImageOutput3D<T>::writeData(
    ScalarField3D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    std::unique_ptr<ScalarField3D<TConv> > transformedField = copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    writeHeader(transformedField->getNx(), transformedField->getNy(), transformedField->getNz());
    vtkOut.writeDataField<TConv>(*transformedField, scalarFieldName, 1);
}

template <typename T>
template <typename TConv>
void VtkAsciiImageOutput3D<T>::writeData(
    MultiScalarField3D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    std::unique_ptr<MultiScalarField3D<TConv> > transformedField =
        copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    writeHeader(transformedField->getBoundingBox());
    vtkOut.writeDataField<TConv>(*transformedField, scalarFieldName, 1);
}

template <typename T>
template <plint n, typename TConv>
void VtkAsciiImageOutput3D<T>::writeData(
    TensorField3D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    std::unique_ptr<TensorField3D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    writeHeader(transformedField->getNx(), transformedField->getNy(), transformedField->getNz());
    vtkOut.writeDataField<TConv>(*transformedField, tensorFieldName, n);
}

template <typename T>
template <plint n, typename TConv>
void VtkAsciiImageOutput3D<T>::writeData(
    MultiTensorField3D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    std::unique_ptr<MultiTensorField3D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    writeHeader(transformedField->getBoundingBox());
    vtkOut.writeDataField<TConv>(*transformedField, tensorFieldName, n);
}

template <typename T>
template <typename TConv>
void VtkAsciiImageOutput3D<T>::writeData(
    MultiNTensorField3D<T> &nTensorField, std::string nTensorFieldName)
{
    MultiNTensorField3D<TConv> *transformedField =
        copyConvert<T, TConv>(nTensorField, nTensorField.getBoundingBox());
    writeHeader(transformedField->getBoundingBox());
    vtkOut.writeDataField<TConv>(*transformedField, nTensorFieldName, transformedField->getNdim());
    delete transformedField;
}

////////// class ParallelVtkImageOutput3D ////////////////////////////////////

template <typename T>
ParallelVtkImageOutput3D<T>::ParallelVtkImageOutput3D(
    std::string fName, plint numEntries_, double deltaX_) :
    fullName(global::directories().getVtkOutDir() + fName + ".vti"),
    vtkOut(fullName),
    deltaX(deltaX_),
    offset(T(), T(), T()),
    headerWritten(false),
    numEntries(numEntries_)
{
    sizeOfEntry = 256;  // Entries are smaller than that. We just make sure to
                        // provide enough space.
    sizeOfFooter = 100;
    nextDataOffset = 0;
    iEntry = 0;
}

template <typename T>
ParallelVtkImageOutput3D<T>::ParallelVtkImageOutput3D(
    std::string fName, plint numEntries_, double deltaX_, Array<double, 3> offset_) :
    fullName(global::directories().getVtkOutDir() + fName + ".vti"),
    vtkOut(fullName),
    deltaX(deltaX_),
    offset(offset_),
    headerWritten(false),
    numEntries(numEntries_)
{
    sizeOfEntry = 256;  // Entries are smaller than that. We just make sure to
                        // provide enough space.
    sizeOfFooter = 100;
    nextDataOffset = 0;
    iEntry = 0;
}

template <typename T>
ParallelVtkImageOutput3D<T>::~ParallelVtkImageOutput3D()
{
    vtkOut.writeFooter2();
}

template <typename T>
void ParallelVtkImageOutput3D<T>::writeHeader(Box3D boundingBox_)
{
    if (headerWritten) {
        PLB_PRECONDITION(boundingBox == boundingBox_);
    } else {
        boundingBox = boundingBox_;
        vtkOut.writeHeader(boundingBox, offset, deltaX);
        vtkOut.startPiece(boundingBox);
        plint fillerSpace = numEntries * sizeOfEntry + sizeOfFooter;
        vtkOut.appendFillerSpace(fillerSpace);
        vtkOut.endPiece();
        vtkOut.writeFooter1();
        headerWritten = true;
    }
}

template <typename T>
template <typename TConv>
void ParallelVtkImageOutput3D<T>::writeData(
    Box3D boundingBox, plint nDim, MultiBlock3D &block, IndexOrdering::OrderingT ordering,
    std::string const &name)
{
    PLB_ASSERT(iEntry < numEntries);

    writeHeader(boundingBox);
    plint declarationPos = vtkOut.getSizeOfHeader() + iEntry * sizeOfEntry;
    vtkOut.declareDataField<TConv>(name, nDim, declarationPos, nextDataOffset);
    plint nextDataSize = block.getBoundingBox().nCells() * nDim * sizeof(TConv);
    nextDataOffset += nextDataSize + sizeof(pluint);  // For the size indicator.
    vtkOut.writeDataField(block, ordering, sizeof(TConv), nextDataSize);
    ++iEntry;
}

template <typename T>
template <typename TConv>
void ParallelVtkImageOutput3D<T>::writeData(
    MultiScalarField3D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    std::unique_ptr<MultiScalarField3D<TConv> > transformedField =
        copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    writeData<TConv>(
        scalarField.getBoundingBox(), 1, *transformedField, IndexOrdering::backward,
        scalarFieldName);
}

template <typename T>
template <plint n, typename TConv>
void ParallelVtkImageOutput3D<T>::writeData(
    MultiTensorField3D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    std::unique_ptr<MultiTensorField3D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    writeData<TConv>(
        transformedField->getBoundingBox(), n, *transformedField, IndexOrdering::backward,
        tensorFieldName);
}

template <typename T>
template <typename TConv>
void ParallelVtkImageOutput3D<T>::writeData(
    MultiNTensorField3D<T> &nTensorField, std::string nTensorFieldName)
{
    MultiNTensorField3D<TConv> *transformedField =
        copyConvert<T, TConv>(nTensorField, nTensorField.getBoundingBox());
    writeData<TConv>(
        transformedField->getBoundingBox(), transformedField->getNdim(), *transformedField,
        IndexOrdering::backward, nTensorFieldName);
    delete transformedField;
}

}  // namespace plb

#endif  // VTK_DATA_OUTPUT_HH
