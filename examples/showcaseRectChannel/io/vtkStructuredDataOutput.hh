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

#ifndef VTK_STRUCTURED_DATA_OUTPUT_HH
#define VTK_STRUCTURED_DATA_OUTPUT_HH

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "core/globalDefs.h"
#include "io/serializerIO.h"
#include "io/vtkDataOutput.hh"
#include "io/vtkStructuredDataOutput.h"
#include "parallelism/mpiManager.h"

namespace plb {

////////// class VtkStructuredWriter3D ////////////////////////////////////////

template <typename T>
void VtkStructuredWriter3D::writeDataField(
    DataSerializer const *serializer, std::string const &name, plint nDim)
{
    if (global::mpi().isMainProcessor()) {
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
    //
    // Those properties are properly handled by the serializer2ostr function, if pluint plint is
    // equal to UInt32. If not, you are on your own.

    bool enforceUint =
        true;  // VTK uses "unsigned" to indicate the size of data, even on a 64-bit machine.
    serializerToBase64Stream(serializer, ostr, enforceUint);

    if (global::mpi().isMainProcessor()) {
        (*ostr) << "\n</DataArray>\n";
    }
}

////////// class VtkStructuredImageOutput2D ////////////////////////////////////

template <typename T>
VtkStructuredImageOutput2D<T>::VtkStructuredImageOutput2D(std::string fName, T deltaX_) :
    fullName(global::directories().getVtkOutDir() + fName + ".vts"),
    vtkOut(fullName),
    deltaX(deltaX_),
    offset(T(), T()),
    headerWritten(false)
{ }

template <typename T>
VtkStructuredImageOutput2D<T>::VtkStructuredImageOutput2D(
    std::string fName, T deltaX_, Array<T, 2> offset_) :
    fullName(global::directories().getVtkOutDir() + fName + ".vts"),
    vtkOut(fullName),
    deltaX(deltaX_),
    offset(offset_),
    headerWritten(false)
{ }

template <typename T>
VtkStructuredImageOutput2D<T>::~VtkStructuredImageOutput2D()
{
    writeFooter();
}

template <typename T>
void VtkStructuredImageOutput2D<T>::writeHeader(plint nx_, plint ny_)
{
    if (headerWritten) {
        PLB_PRECONDITION(nx == nx_);
        PLB_PRECONDITION(ny == ny_);
    } else {
        nx = nx_;
        ny = ny_;
        vtkOut.writeHeader(Box3D(0, nx - 1, 0, ny - 1, 0, 0));
        vtkOut.startPiece(
            Box3D(0, nx - 1, 0, ny - 1, 0, 0), Array<T, 3>(offset[0], offset[1], T()), deltaX);
        headerWritten = true;
    }
}

template <typename T>
void VtkStructuredImageOutput2D<T>::writeFooter()
{
    if (headerWritten) {
        vtkOut.endPiece();
        vtkOut.writeFooter();
        headerWritten = false;
    }
}

template <typename T>
template <typename TConv>
void VtkStructuredImageOutput2D<T>::writeData(
    ScalarField2D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    writeHeader(scalarField.getNx(), scalarField.getNy());
    std::unique_ptr<ScalarField2D<TConv> > transformedField = copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    vtkOut.writeDataField<TConv>(
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        scalarFieldName, 1);
}

template <typename T>
template <typename TConv>
void VtkStructuredImageOutput2D<T>::writeData(
    MultiScalarField2D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    writeHeader(scalarField.getNx(), scalarField.getNy());
    std::unique_ptr<MultiScalarField2D<TConv> > transformedField =
        copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    vtkOut.writeDataField<TConv>(
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        scalarFieldName, 1);
}

template <typename T>
template <plint n, typename TConv>
void VtkStructuredImageOutput2D<T>::writeData(
    TensorField2D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    writeHeader(tensorField.getNx(), tensorField.getNy());
    std::unique_ptr<TensorField2D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    vtkOut.writeDataField<TConv>(
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        tensorFieldName, n);
}

template <typename T>
template <plint n, typename TConv>
void VtkStructuredImageOutput2D<T>::writeData(
    MultiTensorField2D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    writeHeader(tensorField.getNx(), tensorField.getNy());
    std::unique_ptr<MultiTensorField2D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    vtkOut.writeDataField<TConv>(
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        tensorFieldName, n);
}

////////// class VtkStructuredImageOutput3D ////////////////////////////////////

template <typename T>
VtkStructuredImageOutput3D<T>::VtkStructuredImageOutput3D(std::string fName, T deltaX_) :
    fullName(global::directories().getVtkOutDir() + fName + ".vts"),
    vtkOut(fullName),
    deltaX(deltaX_),
    offset(T(), T(), T()),
    headerWritten(false)
{ }

template <typename T>
VtkStructuredImageOutput3D<T>::VtkStructuredImageOutput3D(
    std::string fName, T deltaX_, Array<T, 3> offset_) :
    fullName(global::directories().getVtkOutDir() + fName + ".vts"),
    vtkOut(fullName),
    deltaX(deltaX_),
    offset(offset_),
    headerWritten(false)
{ }

template <typename T>
VtkStructuredImageOutput3D<T>::~VtkStructuredImageOutput3D()
{
    writeFooter();
}

template <typename T>
void VtkStructuredImageOutput3D<T>::writeHeader(plint nx_, plint ny_, plint nz_)
{
    if (headerWritten) {
        PLB_PRECONDITION(nx == nx_);
        PLB_PRECONDITION(ny == ny_);
        PLB_PRECONDITION(nz == nz_);
    } else {
        nx = nx_;
        ny = ny_;
        nz = nz_;
        vtkOut.writeHeader(Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1));
        vtkOut.startPiece(Box3D(0, nx - 1, 0, ny - 1, 0, nz - 1), offset, deltaX);
        headerWritten = true;
    }
}

template <typename T>
void VtkStructuredImageOutput3D<T>::writeFooter()
{
    if (headerWritten) {
        vtkOut.endPiece();
        vtkOut.writeFooter();
        headerWritten = false;
    }
}

template <typename T>
template <typename TConv>
void VtkStructuredImageOutput3D<T>::writeData(
    ScalarField3D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    writeHeader(scalarField.getNx(), scalarField.getNy(), scalarField.getNz());
    std::unique_ptr<ScalarField3D<TConv> > transformedField = copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    vtkOut.writeDataField<TConv>(
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        scalarFieldName, 1);
}

template <typename T>
template <typename TConv>
void VtkStructuredImageOutput3D<T>::writeData(
    MultiScalarField3D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor,
    TConv additiveOffset)
{
    writeHeader(scalarField.getNx(), scalarField.getNy(), scalarField.getNz());
    std::unique_ptr<MultiScalarField3D<TConv> > transformedField =
        copyConvert<T, TConv>(scalarField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    if (!util::isZero(additiveOffset)) {
        addInPlace(*transformedField, additiveOffset);
    }
    vtkOut.writeDataField<TConv>(
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        scalarFieldName, 1);
}

template <typename T>
template <plint n, typename TConv>
void VtkStructuredImageOutput3D<T>::writeData(
    TensorField3D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    writeHeader(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    std::unique_ptr<TensorField3D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    vtkOut.writeDataField<TConv>(
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        tensorFieldName, n);
}

template <typename T>
template <plint n, typename TConv>
void VtkStructuredImageOutput3D<T>::writeData(
    MultiTensorField3D<T, n> &tensorField, std::string tensorFieldName, TConv scalingFactor)
{
    writeHeader(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    std::unique_ptr<MultiTensorField3D<TConv, n> > transformedField =
        copyConvert<T, TConv, n>(tensorField);
    if (!util::isOne(scalingFactor)) {
        multiplyInPlace(*transformedField, scalingFactor);
    }
    vtkOut.writeDataField<TConv>(
        transformedField->getBlockSerializer(
            transformedField->getBoundingBox(), IndexOrdering::backward),
        tensorFieldName, n);
}

}  // namespace plb

#endif  // VTK_STRUCTURED_DATA_OUTPUT_HH
