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

#ifndef VTK_STRUCTURED_DATA_OUTPUT_H
#define VTK_STRUCTURED_DATA_OUTPUT_H

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

class VtkStructuredWriter3D {
public:
    VtkStructuredWriter3D(std::string const &fileName_);
    ~VtkStructuredWriter3D();
    void writeHeader(Box3D domain);
    void startPiece(Box3D domain, const Array<double, 3> &origin, double deltaX);
    void endPiece();
    void writeFooter();
    template <typename T>
    void writeDataField(DataSerializer const *serializer, std::string const &name, plint nDim);

private:
    VtkStructuredWriter3D(VtkStructuredWriter3D const &rhs);
    VtkStructuredWriter3D operator=(VtkStructuredWriter3D const &rhs);

private:
    std::string fileName;
    std::ofstream *ostr;
};

template <typename T>
class VtkStructuredImageOutput2D {
public:
    VtkStructuredImageOutput2D(std::string fName, T deltaX_ = (T)1);
    VtkStructuredImageOutput2D(std::string fName, T deltaX_, Array<T, 2> offset);
    ~VtkStructuredImageOutput2D();
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

private:
    void writeHeader(plint nx_, plint ny_);
    void writeFooter();

private:
    std::string fullName;
    VtkStructuredWriter3D vtkOut;
    T deltaX;
    Array<T, 2> offset;
    bool headerWritten;
    plint nx, ny;
};

template <typename T>
class VtkStructuredImageOutput3D {
public:
    VtkStructuredImageOutput3D(std::string fName, T deltaX_ = (T)1);
    VtkStructuredImageOutput3D(std::string fName, T deltaX_, Array<T, 3> offset);
    ~VtkStructuredImageOutput3D();
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

private:
    void writeHeader(plint nx_, plint ny_, plint nz_);
    void writeFooter();

private:
    std::string fullName;
    VtkStructuredWriter3D vtkOut;
    T deltaX;
    Array<T, 3> offset;
    bool headerWritten;
    plint nx, ny, nz;
};

}  // namespace plb

#endif  // VTK_STRUCTURED_DATA_OUTPUT_H
