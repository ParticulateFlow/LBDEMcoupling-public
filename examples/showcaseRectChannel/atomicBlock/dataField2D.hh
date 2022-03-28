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
 * Scalar, vector and tensor fields for 2D data analysis -- generic implementation.
 */

#ifndef DATA_FIELD_2D_HH
#define DATA_FIELD_2D_HH

#include <algorithm>
#include <cstring>
#include <typeinfo>

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/dataField2D.h"

namespace plb {

/////// Class ScalarField2D //////////////////////////////////

template <typename T>
ScalarField2D<T>::ScalarField2D(plint nx_, plint ny_, T iniVal) :
    AtomicBlock2D(nx_, ny_), dataTransfer(*this)
{
    allocateMemory();
    for (pluint iData = 0; iData < getSize(); ++iData) {
        (*this)[iData] = iniVal;
    }
}

template <typename T>
ScalarField2D<T>::~ScalarField2D()
{
    releaseMemory();
}

template <typename T>
ScalarField2D<T>::ScalarField2D(ScalarField2D<T> const &rhs) :
    AtomicBlock2D(rhs), dataTransfer(*this)
{
    allocateMemory();
    for (pluint iData = 0; iData < getSize(); ++iData) {
        (*this)[iData] = rhs[iData];
    }
}

template <typename T>
ScalarField2D<T> &ScalarField2D<T>::operator=(ScalarField2D<T> const &rhs)
{
    ScalarField2D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T>
void ScalarField2D<T>::swap(ScalarField2D<T> &rhs)
{
    AtomicBlock2D::swap(rhs);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
}

template <typename T>
void ScalarField2D<T>::reset()
{
    for (plint index = 0; index < this->getNx() * this->getNy(); ++index) {
        (*this)[index] = T();
    }
}

template <typename T>
BlockDataTransfer2D &ScalarField2D<T>::getDataTransfer()
{
    return dataTransfer;
}

template <typename T>
BlockDataTransfer2D const &ScalarField2D<T>::getDataTransfer() const
{
    return dataTransfer;
}

template <typename T>
void ScalarField2D<T>::allocateMemory()
{
    rawData = new T[(pluint)this->getNx() * (pluint)this->getNy()];
    field = new T *[(pluint)this->getNx()];
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        field[iX] = rawData + (pluint)iX * (pluint)this->getNy();
    }
}

template <typename T>
void ScalarField2D<T>::releaseMemory()
{
    delete[] field;
    delete[] rawData;
    rawData = 0;
}

////////////////////// Class ScalarFieldDataTransfer2D /////////////////////////

template <typename T>
ScalarFieldDataTransfer2D<T>::ScalarFieldDataTransfer2D(ScalarField2D<T> &field_) : field(field_)
{ }

template <typename T>
plint ScalarFieldDataTransfer2D<T>::staticCellSize() const
{
    return sizeof(T);
}

template <typename T>
void ScalarFieldDataTransfer2D<T>::send(
    Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    PLB_PRECONDITION(contained(domain, field.getBoundingBox()));
    plint cellSize = staticCellSize();
    pluint numBytes = domain.nCells() * cellSize;
    // Avoid dereferencing uninitialized pointer.
    if (numBytes == 0)
        return;
    buffer.resize(numBytes);

    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            memcpy((void *)(&buffer[iData]), (const void *)(&field.get(iX, iY)), sizeof(T));
            iData += sizeof(T);
        }
    }
}

template <typename T>
void ScalarFieldDataTransfer2D<T>::receive(
    Box2D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(contained(domain, field.getBoundingBox()));
    PLB_PRECONDITION(domain.nCells() * staticCellSize() == (plint)buffer.size());

    // Avoid dereferencing uninitialized pointer.
    if (buffer.empty())
        return;
    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            memcpy((void *)(&field.get(iX, iY)), (const void *)(&buffer[iData]), sizeof(T));
            iData += sizeof(T);
        }
    }
}

template <typename T>
void ScalarFieldDataTransfer2D<T>::attribute(
    Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind)
{
    PLB_PRECONDITION(typeid(from) == typeid(ScalarField2D<T> const &));
    PLB_PRECONDITION(contained(toDomain, field.getBoundingBox()));
    ScalarField2D<T> const &fromField = (ScalarField2D<T> const &)from;
    for (plint iX = toDomain.x0; iX <= toDomain.x1; ++iX) {
        for (plint iY = toDomain.y0; iY <= toDomain.y1; ++iY) {
            field.get(iX, iY) = fromField.get(iX + deltaX, iY + deltaY);
        }
    }
}

//////// Class TensorField2D //////////////////////////////////

template <typename T, int nDim>
TensorField2D<T, nDim>::TensorField2D(plint nx_, plint ny_) :
    AtomicBlock2D(nx_, ny_), dataTransfer(*this)
{
    allocateMemory();
    for (plint iData = 0; iData < this->getNx() * this->getNy(); ++iData) {
        for (int iDim = 0; iDim < nDim; ++iDim) {
            (*this)[iData][iDim] = T();
        }
    }
}

template <typename T, int nDim>
TensorField2D<T, nDim>::TensorField2D(plint nx_, plint ny_, Array<T, nDim> const &iniVal) :
    AtomicBlock2D(nx_, ny_), dataTransfer(*this)
{
    allocateMemory();
    for (plint iData = 0; iData < this->getNx() * this->getNy(); ++iData) {
        for (int iDim = 0; iDim < nDim; ++iDim) {
            (*this)[iData][iDim] = iniVal[iDim];
        }
    }
}

template <typename T, int nDim>
TensorField2D<T, nDim>::~TensorField2D()
{
    releaseMemory();
}

template <typename T, int nDim>
TensorField2D<T, nDim>::TensorField2D(TensorField2D<T, nDim> const &rhs) :
    AtomicBlock2D(rhs), dataTransfer(*this)
{
    allocateMemory();
    for (plint iData = 0; iData < this->getNx() * this->getNy(); ++iData) {
        for (int iDim = 0; iDim < nDim; ++iDim) {
            (*this)[iData][iDim] = rhs[iData][iDim];
        }
    }
}

template <typename T, int nDim>
TensorField2D<T, nDim> &TensorField2D<T, nDim>::operator=(TensorField2D<T, nDim> const &rhs)
{
    TensorField2D<T, nDim> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T, int nDim>
void TensorField2D<T, nDim>::swap(TensorField2D<T, nDim> &rhs)
{
    AtomicBlock2D::swap(rhs);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
}

template <typename T, int nDim>
void TensorField2D<T, nDim>::reset()
{
    for (plint index = 0; index < this->getNx() * this->getNy(); ++index) {
        for (int iDim = 0; iDim < nDim; ++iDim) {
            (*this)[index][iDim] = T();
        }
    }
}

template <typename T, int nDim>
BlockDataTransfer2D &TensorField2D<T, nDim>::getDataTransfer()
{
    return dataTransfer;
}

template <typename T, int nDim>
BlockDataTransfer2D const &TensorField2D<T, nDim>::getDataTransfer() const
{
    return dataTransfer;
}

template <typename T, int nDim>
void TensorField2D<T, nDim>::allocateMemory()
{
    rawData = new Array<T, nDim>[(pluint)this->getNx() * (pluint)this->getNy()];
    field = new Array<T, nDim> *[(pluint)this->getNx()];
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        field[iX] = rawData + (pluint)iX * (pluint)this->getNy();
    }
}

template <typename T, int nDim>
void TensorField2D<T, nDim>::releaseMemory()
{
    delete[] field;
    delete[] rawData;
    rawData = 0;
}

////////////////////// Class TensorFieldDataTransfer2D /////////////////////////

template <typename T, int nDim>
TensorFieldDataTransfer2D<T, nDim>::TensorFieldDataTransfer2D(TensorField2D<T, nDim> &field_) :
    field(field_)
{ }

template <typename T, int nDim>
plint TensorFieldDataTransfer2D<T, nDim>::staticCellSize() const
{
    return nDim * sizeof(T);
}

template <typename T, int nDim>
void TensorFieldDataTransfer2D<T, nDim>::send(
    Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    PLB_PRECONDITION(contained(domain, field.getBoundingBox()));
    plint cellSize = staticCellSize();
    pluint numBytes = domain.nCells() * cellSize;
    // Avoid dereferencing uninitialized pointer.
    if (numBytes == 0)
        return;
    buffer.resize(numBytes);

    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            memcpy(
                (void *)(&buffer[iData]), (const void *)(&field.get(iX, iY)[0]), nDim * sizeof(T));
            iData += nDim * sizeof(T);
        }
    }
}

template <typename T, int nDim>
void TensorFieldDataTransfer2D<T, nDim>::receive(
    Box2D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(contained(domain, field.getBoundingBox()));
    PLB_PRECONDITION(domain.nCells() * staticCellSize() == (plint)buffer.size());

    // Avoid dereferencing uninitialized pointer.
    if (buffer.empty())
        return;
    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            memcpy(
                (void *)(&field.get(iX, iY)[0]), (const void *)(&buffer[iData]), nDim * sizeof(T));
            iData += nDim * sizeof(T);
        }
    }
}

template <typename T, int nDim>
void TensorFieldDataTransfer2D<T, nDim>::attribute(
    Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind)
{
    PLB_PRECONDITION(typeid(from) == typeid(TensorField2D<T, nDim> const &));
    PLB_PRECONDITION(contained(toDomain, field.getBoundingBox()));
    TensorField2D<T, nDim> const &fromField = (TensorField2D<T, nDim> const &)from;
    for (plint iX = toDomain.x0; iX <= toDomain.x1; ++iX) {
        for (plint iY = toDomain.y0; iY <= toDomain.y1; ++iY) {
            for (int iDim = 0; iDim < nDim; ++iDim) {
                field.get(iX, iY)[iDim] = fromField.get(iX + deltaX, iY + deltaY)[iDim];
            }
        }
    }
}

//////// Class NTensorField2D //////////////////////////////////

template <typename T>
NTensorField2D<T>::NTensorField2D(plint nx_, plint ny_, plint ndim_) :
    NTensorFieldBase2D<T>(ndim_), AtomicBlock2D(nx_, ny_), dataTransfer(*this)
{
    allocateMemory();
    reset();
}

template <typename T>
NTensorField2D<T>::NTensorField2D(plint nx_, plint ny_, plint ndim_, T const *iniVal) :
    NTensorFieldBase2D<T>(ndim_), AtomicBlock2D(nx_, ny_), dataTransfer(*this)
{
    allocateMemory();
    for (plint iData = 0; iData < this->getNx() * this->getNy() * this->getNdim();
         iData += this->getNdim())
    {
        for (plint iDim = 0; iDim < this->getNdim(); ++iDim) {
            (*this)[iData + iDim] = iniVal[iDim];
        }
    }
}

template <typename T>
NTensorField2D<T>::~NTensorField2D()
{
    releaseMemory();
}

template <typename T>
NTensorField2D<T>::NTensorField2D(NTensorField2D<T> const &rhs) :
    NTensorFieldBase2D<T>(rhs), AtomicBlock2D(rhs), dataTransfer(*this)
{
    allocateMemory();
    for (plint iData = 0; iData < this->getNx() * this->getNy() * this->getNdim(); ++iData) {
        (*this)[iData] = rhs[iData];
    }
}

template <typename T>
NTensorField2D<T> &NTensorField2D<T>::operator=(NTensorField2D<T> const &rhs)
{
    NTensorField2D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T>
void NTensorField2D<T>::swap(NTensorField2D<T> &rhs)
{
    NTensorFieldBase2D<T>::swap(rhs);
    AtomicBlock2D::swap(rhs);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
}

template <typename T>
void NTensorField2D<T>::reset()
{
    for (plint index = 0; index < this->getNx() * this->getNy() * this->getNdim(); ++index) {
        (*this)[index] = T();
    }
}

template <typename T>
BlockDataTransfer2D &NTensorField2D<T>::getDataTransfer()
{
    return dataTransfer;
}

template <typename T>
BlockDataTransfer2D const &NTensorField2D<T>::getDataTransfer() const
{
    return dataTransfer;
}

template <typename T>
void NTensorField2D<T>::allocateMemory()
{
    rawData = new T[(pluint)this->getNx() * (pluint)this->getNy() * (pluint)this->getNdim()];
    field = new T **[(pluint)this->getNx()];
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        field[iX] = new T *[(pluint)this->getNy()];
        for (plint iY = 0; iY < this->getNy(); ++iY) {
            field[iX][iY] =
                rawData
                + (pluint)this->getNdim() * ((pluint)iY + (pluint)this->getNy() * (pluint)iX);
        }
    }
}

template <typename T>
void NTensorField2D<T>::releaseMemory()
{
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        delete[] field[iX];
    }
    delete[] field;
    delete[] rawData;
    rawData = 0;
}

////////////////////// Class NTensorFieldDataTransfer2D /////////////////////////

template <typename T>
NTensorFieldDataTransfer2D<T>::NTensorFieldDataTransfer2D(NTensorField2D<T> &field_) : field(field_)
{ }

template <typename T>
plint NTensorFieldDataTransfer2D<T>::staticCellSize() const
{
    return field.getNdim() * sizeof(T);
}

template <typename T>
void NTensorFieldDataTransfer2D<T>::send(
    Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    PLB_PRECONDITION(contained(domain, field.getBoundingBox()));
    plint cellSize = staticCellSize();
    pluint numBytes = domain.nCells() * cellSize;
    buffer.resize(numBytes);

    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            memcpy((void *)(&buffer[iData]), (const void *)(&field.get(iX, iY)[0]), cellSize);
            iData += cellSize;
        }
    }
}

template <typename T>
void NTensorFieldDataTransfer2D<T>::receive(
    Box2D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(contained(domain, field.getBoundingBox()));
    PLB_PRECONDITION((pluint)domain.nCells() * staticCellSize() == buffer.size());
    plint cellSize = staticCellSize();

    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            memcpy((void *)(&field.get(iX, iY)[0]), (const void *)(&buffer[iData]), cellSize);
            iData += cellSize;
        }
    }
}

template <typename T>
void NTensorFieldDataTransfer2D<T>::attribute(
    Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind)
{
    PLB_PRECONDITION(typeid(from) == typeid(NTensorField2D<T> const &));
    PLB_PRECONDITION(contained(toDomain, field.getBoundingBox()));
    NTensorField2D<T> const &fromField = (NTensorField2D<T> const &)from;
    for (plint iX = toDomain.x0; iX <= toDomain.x1; ++iX) {
        for (plint iY = toDomain.y0; iY <= toDomain.y1; ++iY) {
            for (int iDim = 0; iDim < field.getNdim(); ++iDim) {
                field.get(iX, iY)[iDim] = fromField.get(iX + deltaX, iY + deltaY)[iDim];
            }
        }
    }
}

}  // namespace plb

#endif  // DATA_FIELD_2D_HH
