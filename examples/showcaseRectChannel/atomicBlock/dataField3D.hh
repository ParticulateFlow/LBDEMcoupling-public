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
 * Scalar, vector and tensor fields for 3D data analysis -- generic implementation.
 */
#ifndef DATA_FIELD_3D_HH
#define DATA_FIELD_3D_HH

#include <algorithm>
#include <cstring>
#include <typeinfo>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/plbTimer.h"

namespace plb {

/////// Class ScalarField3D //////////////////////////////////

template <typename T>
ScalarField3D<T>::ScalarField3D(plint nx_, plint ny_, plint nz_, T iniVal) :
    AtomicBlock3D(nx_, ny_, nz_, new ScalarFieldDataTransfer3D<T>()), ownsMemory(true)
{
    allocateMemory();
    for (pluint iData = 0; iData < getSize(); ++iData) {
        (*this)[iData] = iniVal;
    }
    global::plbCounter("MEMORY_SCALAR").increment(allocatedMemory());
}

template <typename T>
ScalarField3D<T>::~ScalarField3D()
{
    global::plbCounter("MEMORY_SCALAR").increment(-allocatedMemory());
    releaseMemory();
}

template <typename T>
ScalarField3D<T>::ScalarField3D(ScalarField3D<T> const &rhs) : AtomicBlock3D(rhs), ownsMemory(true)
{
    allocateMemory();
    for (pluint iData = 0; iData < getSize(); ++iData) {
        (*this)[iData] = rhs[iData];
    }
    global::plbCounter("MEMORY_SCALAR").increment(allocatedMemory());
}

template <typename T>
ScalarField3D<T>::ScalarField3D(NTensorField3D<T> &rhs) :
    AtomicBlock3D(rhs, new ScalarFieldDataTransfer3D<T>()), ownsMemory(false)
{
    rawData = rhs.rawData;
    allocateMemory();
    global::plbCounter("MEMORY_SCALAR").increment(allocatedMemory());
}

template <typename T>
ScalarField3D<T> &ScalarField3D<T>::operator=(ScalarField3D<T> const &rhs)
{
    ScalarField3D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T>
void ScalarField3D<T>::swap(ScalarField3D<T> &rhs)
{
    global::plbCounter("MEMORY_SCALAR").increment(-allocatedMemory());
    AtomicBlock3D::swap(rhs);
    std::swap(ownsMemory, rhs.ownsMemory);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
    global::plbCounter("MEMORY_SCALAR").increment(allocatedMemory());
}

template <typename T>
void ScalarField3D<T>::reset()
{
    for (plint index = 0; index < this->getNx() * this->getNy() * this->getNz(); ++index) {
        (*this)[index] = T();
    }
}

template <typename T>
void ScalarField3D<T>::allocateMemory()
{
    if (ownsMemory) {
        rawData = new T[(pluint)this->getNx() * (pluint)this->getNy() * (pluint)this->getNz()];
    }
    field = new T **[(pluint)this->getNx()];
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        field[iX] = new T *[(pluint)this->getNy()];
        for (plint iY = 0; iY < this->getNy(); ++iY) {
            field[iX][iY] =
                rawData + (pluint)this->getNz() * ((pluint)iY + (pluint)this->getNy() * (pluint)iX);
        }
    }
}

template <typename T>
void ScalarField3D<T>::releaseMemory()
{
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        delete[] field[iX];
    }
    delete[] field;
    if (ownsMemory) {
        delete[] rawData;
        rawData = 0;
    }
}

template <typename T>
plint ScalarField3D<T>::allocatedMemory() const
{
    if (ownsMemory) {
        return this->getNx() * this->getNy() * this->getNz() * sizeof(T);
    } else {
        return 0;
    }
}

////////////////////// Class ScalarFieldDataTransfer3D /////////////////////////

template <typename T>
ScalarFieldDataTransfer3D<T>::ScalarFieldDataTransfer3D() : field(0), constField(0)
{ }

template <typename T>
void ScalarFieldDataTransfer3D<T>::setBlock(AtomicBlock3D &block)
{
    field = dynamic_cast<ScalarField3D<T> *>(&block);
    PLB_ASSERT(field);
    constField = field;
}

template <typename T>
void ScalarFieldDataTransfer3D<T>::setConstBlock(AtomicBlock3D const &block)
{
    constField = dynamic_cast<ScalarField3D<T> const *>(&block);
    PLB_ASSERT(constField);
}

template <typename T>
ScalarFieldDataTransfer3D<T> *ScalarFieldDataTransfer3D<T>::clone() const
{
    return new ScalarFieldDataTransfer3D<T>(*this);
}

template <typename T>
plint ScalarFieldDataTransfer3D<T>::staticCellSize() const
{
    return sizeof(T);
}

template <typename T>
void ScalarFieldDataTransfer3D<T>::send(
    Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    PLB_PRECONDITION(constField);
    PLB_PRECONDITION(contained(domain, constField->getBoundingBox()));
    plint cellSize = staticCellSize();
    pluint numBytes = domain.nCells() * cellSize;
    // Avoid dereferencing uninitialized pointer.
    if (numBytes == 0)
        return;
    buffer.resize(numBytes);

    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                memcpy(
                    (void *)(&buffer[iData]), (const void *)(&constField->get(iX, iY, iZ)),
                    sizeof(T));
                iData += sizeof(T);
            }
        }
    }
}

template <typename T>
void ScalarFieldDataTransfer3D<T>::receive(
    Box3D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(field);
    PLB_PRECONDITION(contained(domain, field->getBoundingBox()));
    PLB_PRECONDITION(domain.nCells() * staticCellSize() == (plint)buffer.size());

    // Avoid dereferencing uninitialized pointer.
    if (buffer.empty())
        return;
    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                memcpy(
                    (void *)(&field->get(iX, iY, iZ)), (const void *)(&buffer[iData]), sizeof(T));
                iData += sizeof(T);
            }
        }
    }
}

template <typename T>
void ScalarFieldDataTransfer3D<T>::attribute(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
    modif::ModifT kind)
{
    PLB_PRECONDITION(typeid(from) == typeid(ScalarField3D<T> const &));
    PLB_PRECONDITION(contained(toDomain, field->getBoundingBox()));
    ScalarField3D<T> const &fromField = (ScalarField3D<T> const &)from;
    for (plint iX = toDomain.x0; iX <= toDomain.x1; ++iX) {
        for (plint iY = toDomain.y0; iY <= toDomain.y1; ++iY) {
            for (plint iZ = toDomain.z0; iZ <= toDomain.z1; ++iZ) {
                field->get(iX, iY, iZ) = fromField.get(iX + deltaX, iY + deltaY, iZ + deltaZ);
            }
        }
    }
}

//////// Class TensorField3D //////////////////////////////////

template <typename T, int nDim>
TensorField3D<T, nDim>::TensorField3D(plint nx_, plint ny_, plint nz_) :
    AtomicBlock3D(nx_, ny_, nz_, new TensorFieldDataTransfer3D<T, nDim>()), ownsMemory(true)
{
    allocateMemory();
    for (plint iData = 0; iData < this->getNx() * this->getNy() * this->getNz(); ++iData) {
        for (int iDim = 0; iDim < nDim; ++iDim) {
            (*this)[iData][iDim] = T();
        }
    }
    global::plbCounter("MEMORY_TENSOR").increment(allocatedMemory());
}

template <typename T, int nDim>
TensorField3D<T, nDim>::TensorField3D(
    plint nx_, plint ny_, plint nz_, Array<T, nDim> const &iniVal) :
    AtomicBlock3D(nx_, ny_, nz_, new TensorFieldDataTransfer3D<T, nDim>()), ownsMemory(true)
{
    allocateMemory();
    for (plint iData = 0; iData < this->getNx() * this->getNy() * this->getNz(); ++iData) {
        for (int iDim = 0; iDim < nDim; ++iDim) {
            (*this)[iData][iDim] = iniVal[iDim];
        }
    }
    global::plbCounter("MEMORY_TENSOR").increment(allocatedMemory());
}

template <typename T, int nDim>
TensorField3D<T, nDim>::~TensorField3D()
{
    global::plbCounter("MEMORY_TENSOR").increment(-allocatedMemory());
    releaseMemory();
}

template <typename T, int nDim>
TensorField3D<T, nDim>::TensorField3D(TensorField3D<T, nDim> const &rhs) :
    AtomicBlock3D(rhs), ownsMemory(true)
{
    allocateMemory();
    for (plint iData = 0; iData < this->getNx() * this->getNy() * this->getNz(); ++iData) {
        for (int iDim = 0; iDim < nDim; ++iDim) {
            (*this)[iData][iDim] = rhs[iData][iDim];
        }
    }
    global::plbCounter("MEMORY_TENSOR").increment(allocatedMemory());
}

template <typename T, int nDim>
TensorField3D<T, nDim>::TensorField3D(NTensorField3D<T> &rhs) :
    AtomicBlock3D(rhs, new TensorFieldDataTransfer3D<T, nDim>()), ownsMemory(false)
{
    rawData = reinterpret_cast<Array<T, nDim> *>(rhs.rawData);
    allocateMemory();
    global::plbCounter("MEMORY_TENSOR").increment(allocatedMemory());
}

template <typename T, int nDim>
TensorField3D<T, nDim> &TensorField3D<T, nDim>::operator=(TensorField3D<T, nDim> const &rhs)
{
    TensorField3D<T, nDim> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T, int nDim>
void TensorField3D<T, nDim>::swap(TensorField3D<T, nDim> &rhs)
{
    global::plbCounter("MEMORY_TENSOR").increment(-allocatedMemory());
    std::swap(ownsMemory, rhs.ownsMemory);
    AtomicBlock3D::swap(rhs);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
    global::plbCounter("MEMORY_TENSOR").increment(allocatedMemory());
}

template <typename T, int nDim>
void TensorField3D<T, nDim>::reset()
{
    for (plint index = 0; index < this->getNx() * this->getNy() * this->getNz(); ++index) {
        for (int iDim = 0; iDim < nDim; ++iDim) {
            (*this)[index][iDim] = T();
        }
    }
}

template <typename T, int nDim>
void TensorField3D<T, nDim>::allocateMemory()
{
    if (ownsMemory) {
        rawData = new Array<
            T, nDim>[(pluint)this->getNx() * (pluint)this->getNy() * (pluint)this->getNz()];
    }
    field = new Array<T, nDim> **[(pluint)this->getNx()];
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        field[iX] = new Array<T, nDim> *[(pluint)this->getNy()];
        for (plint iY = 0; iY < this->getNy(); ++iY) {
            field[iX][iY] =
                rawData + (pluint)this->getNz() * ((pluint)iY + (pluint)this->getNy() * (pluint)iX);
        }
    }
}

template <typename T, int nDim>
plint TensorField3D<T, nDim>::allocatedMemory() const
{
    if (ownsMemory) {
        return this->getNx() * this->getNy() * this->getNz() * sizeof(T) * nDim;
    } else {
        return 0;
    }
}

template <typename T, int nDim>
void TensorField3D<T, nDim>::releaseMemory()
{
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        delete[] field[iX];
    }
    delete[] field;
    if (ownsMemory) {
        delete[] rawData;
        rawData = 0;
    }
}

////////////////////// Class TensorFieldDataTransfer3D /////////////////////////

template <typename T, int nDim>
TensorFieldDataTransfer3D<T, nDim>::TensorFieldDataTransfer3D() : field(0), constField(0)
{ }

template <typename T, int nDim>
void TensorFieldDataTransfer3D<T, nDim>::setBlock(AtomicBlock3D &block)
{
    field = dynamic_cast<TensorField3D<T, nDim> *>(&block);
    PLB_ASSERT(field);
    constField = field;
}

template <typename T, int nDim>
void TensorFieldDataTransfer3D<T, nDim>::setConstBlock(AtomicBlock3D const &block)
{
    constField = dynamic_cast<TensorField3D<T, nDim> const *>(&block);
    PLB_ASSERT(constField);
}

template <typename T, int nDim>
TensorFieldDataTransfer3D<T, nDim> *TensorFieldDataTransfer3D<T, nDim>::clone() const
{
    return new TensorFieldDataTransfer3D<T, nDim>(*this);
}

template <typename T, int nDim>
plint TensorFieldDataTransfer3D<T, nDim>::staticCellSize() const
{
    return nDim * sizeof(T);
}

template <typename T, int nDim>
void TensorFieldDataTransfer3D<T, nDim>::send(
    Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    PLB_PRECONDITION(constField);
    PLB_PRECONDITION(contained(domain, constField->getBoundingBox()));
    plint cellSize = staticCellSize();
    pluint numBytes = domain.nCells() * cellSize;
    // Avoid dereferencing uninitialized pointer.
    if (numBytes == 0)
        return;
    buffer.resize(numBytes);

    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                memcpy(
                    (void *)(&buffer[iData]), (const void *)(&constField->get(iX, iY, iZ)[0]),
                    nDim * sizeof(T));
                iData += nDim * sizeof(T);
            }
        }
    }
}

template <typename T, int nDim>
void TensorFieldDataTransfer3D<T, nDim>::receive(
    Box3D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(field);
    PLB_PRECONDITION(contained(domain, field->getBoundingBox()));
    PLB_PRECONDITION(domain.nCells() * staticCellSize() == (plint)buffer.size());

    // Avoid dereferencing uninitialized pointer.
    if (buffer.empty())
        return;
    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                memcpy(
                    (void *)(&field->get(iX, iY, iZ)[0]), (const void *)(&buffer[iData]),
                    nDim * sizeof(T));
                iData += nDim * sizeof(T);
            }
        }
    }
}

template <typename T, int nDim>
void TensorFieldDataTransfer3D<T, nDim>::attribute(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
    modif::ModifT kind)
{
    PLB_PRECONDITION(typeid(from) == typeid(TensorField3D<T, nDim> const &));
    PLB_PRECONDITION(contained(toDomain, field->getBoundingBox()));
    TensorField3D<T, nDim> const &fromField = (TensorField3D<T, nDim> const &)from;
    for (plint iX = toDomain.x0; iX <= toDomain.x1; ++iX) {
        for (plint iY = toDomain.y0; iY <= toDomain.y1; ++iY) {
            for (plint iZ = toDomain.z0; iZ <= toDomain.z1; ++iZ) {
                for (int iDim = 0; iDim < nDim; ++iDim) {
                    field->get(iX, iY, iZ)[iDim] =
                        fromField.get(iX + deltaX, iY + deltaY, iZ + deltaZ)[iDim];
                }
            }
        }
    }
}

//////// Class NTensorField3D //////////////////////////////////

template <typename T>
NTensorField3D<T>::NTensorField3D(plint nx_, plint ny_, plint nz_, plint ndim_) :
    NTensorFieldBase3D<T>(ndim_),
    AtomicBlock3D(nx_, ny_, nz_, new NTensorFieldDataTransfer3D<T>()),
    ownsMemory(true)
{
    allocateMemory();
    reset();
    global::plbCounter("MEMORY_NTENSOR").increment(allocatedMemory());
}

template <typename T>
NTensorField3D<T>::NTensorField3D(plint nx_, plint ny_, plint nz_, plint ndim_, T const *iniVal) :
    NTensorFieldBase3D<T>(ndim_),
    AtomicBlock3D(nx_, ny_, nz_, new NTensorFieldDataTransfer3D<T>()),
    ownsMemory(true)
{
    allocateMemory();
    for (plint iData = 0; iData < this->getNx() * this->getNy() * this->getNz() * this->getNdim();
         iData += this->getNdim())
    {
        for (plint iDim = 0; iDim < this->getNdim(); ++iDim) {
            (*this)[iData + iDim] = iniVal[iDim];
        }
    }
    global::plbCounter("MEMORY_NTENSOR").increment(allocatedMemory());
}

template <typename T>
NTensorField3D<T>::NTensorField3D(ScalarField3D<T> &rhs) :
    NTensorFieldBase3D<T>(1),
    AtomicBlock3D(rhs, new NTensorFieldDataTransfer3D<T>()),
    ownsMemory(false)
{
    rawData = rhs.rawData;
    allocateMemory();
    global::plbCounter("MEMORY_NTENSOR").increment(allocatedMemory());
}

template <typename T>
template <int nDim>
NTensorField3D<T>::NTensorField3D(TensorField3D<T, nDim> &rhs) :
    NTensorFieldBase3D<T>(nDim),
    AtomicBlock3D(rhs, new NTensorFieldDataTransfer3D<T>()),
    ownsMemory(false)
{
    rawData = reinterpret_cast<T *>(rhs.rawData);
    allocateMemory();
    global::plbCounter("MEMORY_NTENSOR").increment(allocatedMemory());
}

template <typename T>
NTensorField3D<T>::~NTensorField3D()
{
    global::plbCounter("MEMORY_NTENSOR").increment(-allocatedMemory());
    releaseMemory();
}

template <typename T>
NTensorField3D<T>::NTensorField3D(NTensorField3D<T> const &rhs) :
    NTensorFieldBase3D<T>(rhs), AtomicBlock3D(rhs), ownsMemory(true)
{
    allocateMemory();
    for (plint iData = 0; iData < this->getNx() * this->getNy() * this->getNz() * this->getNdim();
         ++iData)
    {
        (*this)[iData] = rhs[iData];
    }
    global::plbCounter("MEMORY_NTENSOR").increment(allocatedMemory());
}

template <typename T>
NTensorField3D<T> &NTensorField3D<T>::operator=(NTensorField3D<T> const &rhs)
{
    NTensorField3D<T> tmp(rhs);
    swap(tmp);
    return *this;
}

template <typename T>
void NTensorField3D<T>::swap(NTensorField3D<T> &rhs)
{
    global::plbCounter("MEMORY_NTENSOR").increment(-allocatedMemory());
    NTensorFieldBase3D<T>::swap(rhs);
    AtomicBlock3D::swap(rhs);
    std::swap(ownsMemory, rhs.ownsMemory);
    std::swap(rawData, rhs.rawData);
    std::swap(field, rhs.field);
    global::plbCounter("MEMORY_NTENSOR").increment(allocatedMemory());
}

template <typename T>
void NTensorField3D<T>::reset()
{
    for (plint index = 0; index < this->getNx() * this->getNy() * this->getNz() * this->getNdim();
         ++index)
    {
        (*this)[index] = T();
    }
}

template <typename T>
void NTensorField3D<T>::allocateMemory()
{
    if (ownsMemory) {
        rawData = new T
            [(pluint)this->getNx() * (pluint)this->getNy() * (pluint)this->getNz()
             * (pluint)this->getNdim()];
    }
    field = new T ***[(pluint)this->getNx()];
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        field[iX] = new T **[(pluint)this->getNy()];
        for (plint iY = 0; iY < this->getNy(); ++iY) {
            field[iX][iY] = new T *[(plint)this->getNz()];
            for (plint iZ = 0; iZ < this->getNz(); ++iZ) {
                field[iX][iY][iZ] =
                    rawData
                    + (pluint)this->getNdim()
                          * ((pluint)iZ
                             + (pluint)this->getNz()
                                   * ((pluint)iY + (pluint)this->getNy() * (pluint)iX));
            }
        }
    }
}

template <typename T>
plint NTensorField3D<T>::allocatedMemory() const
{
    if (ownsMemory) {
        return this->getNx() * this->getNy() * this->getNz() * sizeof(T) * this->getNdim();
    } else {
        return 0;
    }
}

template <typename T>
void NTensorField3D<T>::releaseMemory()
{
    for (plint iX = 0; iX < this->getNx(); ++iX) {
        for (plint iY = 0; iY < this->getNy(); ++iY) {
            delete[] field[iX][iY];
        }
        delete[] field[iX];
    }
    delete[] field;
    if (ownsMemory) {
        delete[] rawData;
        rawData = 0;
    }
}

////////////////////// Class NTensorFieldDataTransfer3D /////////////////////////

template <typename T>
NTensorFieldDataTransfer3D<T>::NTensorFieldDataTransfer3D() : field(0), constField(0)
{ }

template <typename T>
void NTensorFieldDataTransfer3D<T>::setBlock(AtomicBlock3D &block)
{
    field = dynamic_cast<NTensorField3D<T> *>(&block);
    PLB_ASSERT(field);
    constField = field;
}

template <typename T>
void NTensorFieldDataTransfer3D<T>::setConstBlock(AtomicBlock3D const &block)
{
    constField = dynamic_cast<NTensorField3D<T> const *>(&block);
    PLB_ASSERT(constField);
}

template <typename T>
NTensorFieldDataTransfer3D<T> *NTensorFieldDataTransfer3D<T>::clone() const
{
    return new NTensorFieldDataTransfer3D<T>(*this);
}

template <typename T>
plint NTensorFieldDataTransfer3D<T>::staticCellSize() const
{
    return constField->getNdim() * sizeof(T);
}

template <typename T>
void NTensorFieldDataTransfer3D<T>::send(
    Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const
{
    PLB_PRECONDITION(constField);
    PLB_PRECONDITION(contained(domain, constField->getBoundingBox()));
    plint cellSize = staticCellSize();
    pluint numBytes = domain.nCells() * cellSize;
    buffer.resize(numBytes);

    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                memcpy(
                    (void *)(&buffer[iData]), (const void *)(&constField->get(iX, iY, iZ)[0]),
                    cellSize);
                iData += cellSize;
            }
        }
    }
}

template <typename T>
void NTensorFieldDataTransfer3D<T>::receive(
    Box3D domain, std::vector<char> const &buffer, modif::ModifT kind)
{
    PLB_PRECONDITION(field);
    PLB_PRECONDITION(contained(domain, field->getBoundingBox()));
    PLB_PRECONDITION((pluint)domain.nCells() * staticCellSize() == buffer.size());
    plint cellSize = staticCellSize();

    plint iData = 0;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                memcpy(
                    (void *)(&field->get(iX, iY, iZ)[0]), (const void *)(&buffer[iData]), cellSize);
                iData += cellSize;
            }
        }
    }
}

template <typename T>
void NTensorFieldDataTransfer3D<T>::attribute(
    Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
    modif::ModifT kind)
{
    PLB_PRECONDITION(typeid(from) == typeid(NTensorField3D<T> const &));
    PLB_PRECONDITION(contained(toDomain, field->getBoundingBox()));
    NTensorField3D<T> const &fromField = (NTensorField3D<T> const &)from;
    for (plint iX = toDomain.x0; iX <= toDomain.x1; ++iX) {
        for (plint iY = toDomain.y0; iY <= toDomain.y1; ++iY) {
            for (plint iZ = toDomain.z0; iZ <= toDomain.z1; ++iZ) {
                for (int iDim = 0; iDim < field->getNdim(); ++iDim) {
                    field->get(iX, iY, iZ)[iDim] =
                        fromField.get(iX + deltaX, iY + deltaY, iZ + deltaZ)[iDim];
                }
            }
        }
    }
}

}  // namespace plb

#endif  // DATA_FIELD_3D_HH
