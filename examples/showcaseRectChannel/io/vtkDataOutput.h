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

#ifndef VTK_DATA_OUTPUT_H
#define VTK_DATA_OUTPUT_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataField3D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "core/serializer.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

class VtkDataWriter3D {
public:
    VtkDataWriter3D(
        std::string const &fileName_, bool pointData_ = true, bool mainProcOnly_ = true);
    ~VtkDataWriter3D();
    void writeHeader(Box3D domain, Array<double, 3> origin, double deltaX);
    void startPiece(Box3D domain);
    void endPiece();
    void writeFooter();
    template <typename T>
    void writeDataField(DataSerializer const *serializer, std::string const &name, plint nDim);

private:
    VtkDataWriter3D(VtkDataWriter3D const &rhs);
    VtkDataWriter3D operator=(VtkDataWriter3D const &rhs);

private:
    std::string fileName;
    std::ofstream *ostr;
    bool pointData;
    bool mainProcOnly;
};

class VtkAsciiDataWriter3D {
public:
    VtkAsciiDataWriter3D(
        std::string const &fileName_, bool pointData_ = true, bool mainProcOnly_ = true);
    ~VtkAsciiDataWriter3D();
    void writeHeader(Box3D domain, Array<double, 3> origin, double deltaX);
    void startPiece(Box3D domain, std::string const &scalars);
    void endPiece();
    void writeFooter();
    template <typename T>
    void writeDataField(Block3D const &block, std::string const &name, plint nDim);

private:
    VtkAsciiDataWriter3D(VtkAsciiDataWriter3D const &rhs);
    VtkAsciiDataWriter3D operator=(VtkAsciiDataWriter3D const &rhs);

private:
    std::string fileName;
    std::ofstream *ostr;
    bool pointData;
    bool mainProcOnly;
};

class ParallelVtkDataWriter3D {
public:
    ParallelVtkDataWriter3D(std::string const &fileName_);
    ~ParallelVtkDataWriter3D();
    void writeHeader(Box3D domain, Array<double, 3> origin, double deltaX);
    void startPiece(Box3D domain);
    void endPiece();
    void writeFooter1();
    void writeFooter2();
    template <typename T>
    void declareDataField(
        std::string const &name, plint nDim, plint declarationPos, plint dataOffset);
    void writeDataField(
        MultiBlock3D &block, IndexOrdering::OrderingT ordering, plint sizeOfScalar, plint dataSize);
    void appendFillerSpace(plint fillerSpace);
    plint getSizeOfHeader() const
    {
        return sizeOfHeader;
    }

private:
    ParallelVtkDataWriter3D(ParallelVtkDataWriter3D const &rhs);
    ParallelVtkDataWriter3D operator=(ParallelVtkDataWriter3D const &rhs);

private:
    plint sizeOfHeader;
    std::string fileName;
    std::ofstream *ostr;
};

template <typename T>
class VtkImageOutput2D {
public:
    VtkImageOutput2D(std::string fName, double deltaX_ = 1.);
    VtkImageOutput2D(std::string fName, double deltaX_, Array<double, 2> offset);
    ~VtkImageOutput2D();
    template <typename TConv>
    void writeData(
        plint nx, plint ny, plint nDim, DataSerializer const *serializer, std::string const &name);
    template <typename TConv>
    void writeData(
        Box2D boundingBox_, plint nDim, DataSerializer const *serializer, std::string const &name);
    template <typename TConv>
    void writeData(
        ScalarField2D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor = (TConv)1,
        TConv additiveOffset = (TConv)0);
    template <typename TConv>
    void writeData(
        MultiScalarField2D<T> &scalarField, std::string scalarFieldName,
        TConv scalingFactor = (TConv)1, TConv additiveOffset = (TConv)0);
    template <plint n, typename TConv>
    void writeData(
        TensorField2D<T, n> &tensorField, std::string tensorFieldName,
        TConv scalingFactor = (TConv)1);
    template <plint n, typename TConv>
    void writeData(
        MultiTensorField2D<T, n> &tensorField, std::string tensorFieldName,
        TConv scalingFactor = (TConv)1);
    template <typename TConv>
    void writeData(MultiNTensorField2D<T> &nTensorField, std::string nTensorFieldName);

private:
    void writeHeader(plint nx_, plint ny_);
    void writeHeader(Box2D boundingBox_);
    void writeFooter();

private:
    std::string fullName;
    VtkDataWriter3D vtkOut;
    double deltaX;
    Array<T, 2> offset;
    bool headerWritten;
    plint nx, ny;
    Box3D boundingBox;
};

template <typename T>
class VtkImageOutput3D {
public:
    VtkImageOutput3D(std::string fName, double deltaX_ = 1.);
    VtkImageOutput3D(std::string fName, double deltaX_, Array<double, 3> offset);
    ~VtkImageOutput3D();
    template <typename TConv>
    void writeData(
        plint nx, plint ny, plint nz, plint nDim, DataSerializer const *serializer,
        std::string const &name);
    template <typename TConv>
    void writeData(
        Box3D boundingBox_, plint nDim, DataSerializer const *serializer, std::string const &name);
    template <typename TConv>
    void writeData(
        ScalarField3D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor = (TConv)1,
        TConv additiveOffset = (TConv)0);
    template <typename TConv>
    void writeData(
        MultiScalarField3D<T> &scalarField, std::string scalarFieldName,
        TConv scalingFactor = (TConv)1, TConv additiveOffset = (TConv)0);
    template <plint n, typename TConv>
    void writeData(
        TensorField3D<T, n> &tensorField, std::string tensorFieldName,
        TConv scalingFactor = (TConv)1);
    template <plint n, typename TConv>
    void writeData(
        MultiTensorField3D<T, n> &tensorField, std::string tensorFieldName,
        TConv scalingFactor = (TConv)1);
    template <typename TConv>
    void writeData(MultiNTensorField3D<T> &nTensorField, std::string nTensorFieldName);

private:
    void writeHeader(plint nx_, plint ny_, plint nz_);
    void writeHeader(Box3D boundingBox_);
    void writeFooter();

private:
    std::string fullName;
    VtkDataWriter3D vtkOut;
    double deltaX;
    Array<T, 3> offset;
    bool headerWritten;
    Box3D boundingBox;
};

template <typename T>
class VtkAsciiImageOutput3D {
public:
    VtkAsciiImageOutput3D(std::string fName, std::string const &scalars_, double deltaX_ = 1.);
    VtkAsciiImageOutput3D(
        std::string fName, std::string const &scalars_, double deltaX_, Array<double, 3> offset);
    ~VtkAsciiImageOutput3D();
    template <typename TConv>
    void writeData(
        ScalarField3D<T> &scalarField, std::string scalarFieldName, TConv scalingFactor = (TConv)1,
        TConv additiveOffset = (TConv)0);
    template <typename TConv>
    void writeData(
        MultiScalarField3D<T> &scalarField, std::string scalarFieldName,
        TConv scalingFactor = (TConv)1, TConv additiveOffset = (TConv)0);
    template <plint n, typename TConv>
    void writeData(
        TensorField3D<T, n> &tensorField, std::string tensorFieldName,
        TConv scalingFactor = (TConv)1);
    template <plint n, typename TConv>
    void writeData(
        MultiTensorField3D<T, n> &tensorField, std::string tensorFieldName,
        TConv scalingFactor = (TConv)1);
    template <typename TConv>
    void writeData(MultiNTensorField3D<T> &nTensorField, std::string nTensorFieldName);

private:
    void writeHeader(plint nx_, plint ny_, plint nz_);
    void writeHeader(Box3D boundingBox_);
    void writeFooter();

private:
    std::string fullName;
    VtkAsciiDataWriter3D vtkOut;
    std::string scalars;
    double deltaX;
    Array<T, 3> offset;
    bool headerWritten;
    Box3D boundingBox;
};

template <typename T>
class ParallelVtkImageOutput3D {
public:
    ParallelVtkImageOutput3D(std::string fName, plint numEntries, double deltaX_ = 1.);
    ParallelVtkImageOutput3D(
        std::string fName, plint numEntries, double deltaX_, Array<double, 3> offset);
    ~ParallelVtkImageOutput3D();
    template <typename TConv>
    void writeData(
        Box3D boundingBox_, plint nDim, MultiBlock3D &block, IndexOrdering::OrderingT ordering,
        std::string const &name);
    template <typename TConv>
    void writeData(
        MultiScalarField3D<T> &scalarField, std::string scalarFieldName,
        TConv scalingFactor = (TConv)1, TConv additiveOffset = (TConv)0);
    template <plint n, typename TConv>
    void writeData(
        MultiTensorField3D<T, n> &tensorField, std::string tensorFieldName,
        TConv scalingFactor = (TConv)1);
    template <typename TConv>
    void writeData(MultiNTensorField3D<T> &nTensorField, std::string nTensorFieldName);

private:
    void writeHeader(Box3D boundingBox_);

private:
    std::string fullName;
    ParallelVtkDataWriter3D vtkOut;
    double deltaX;
    Array<T, 3> offset;
    bool headerWritten;
    Box3D boundingBox;
    plint numEntries, sizeOfEntry, sizeOfFooter, nextDataOffset, iEntry;
};

}  // namespace plb

#endif  // VTK_DATA_OUTPUT_H
