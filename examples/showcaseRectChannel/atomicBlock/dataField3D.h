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
 * Serial implementation of scalar, vector and tensor fields for 3D data analysis.
 * -- header file
 */

#ifndef DATA_FIELD_3D_H
#define DATA_FIELD_3D_H

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/dataField2D.h"
#include "core/dataFieldBase2D.h"
#include "core/dataFieldBase3D.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"

namespace plb {

template <typename T>
class ScalarField3D;

template <typename T>
class ScalarFieldDataTransfer3D : public BlockDataTransfer3D {
public:
    ScalarFieldDataTransfer3D();
    virtual void setBlock(AtomicBlock3D &block);
    virtual void setConstBlock(AtomicBlock3D const &block);
    virtual ScalarFieldDataTransfer3D<T> *clone() const;
    virtual plint staticCellSize() const;
    /// Send data from the block into a byte-stream.
    virtual void send(Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the block.
    virtual void receive(Box3D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset)
    {
        receive(domain, buffer, kind);
    }
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    {
        receive(domain, buffer, kind);
    }
    /// Attribute data between two blocks.
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind);
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind, Dot3D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
    }

private:
    ScalarField3D<T> *field;
    ScalarField3D<T> const *constField;
};

template <typename T>
class ScalarField3D : public ScalarFieldBase3D<T>, public AtomicBlock3D {
public:
    ScalarField3D(plint nx_, plint ny_, plint nz_, T iniVal = T());
    ~ScalarField3D();
    ScalarField3D(ScalarField3D<T> const &rhs);
    ScalarField3D(NTensorField3D<T> &rhs);
    ScalarField3D<T> &operator=(ScalarField3D<T> const &rhs);
    void swap(ScalarField3D<T> &rhs);

public:
    virtual void reset();
    virtual pluint getSize() const
    {
        return (pluint)this->getNx() * (pluint)this->getNy() * (pluint)this->getNz();
    }
    virtual T &get(plint iX, plint iY, plint iZ)
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
        return field[iX][iY][iZ];
    }
    virtual T const &get(plint iX, plint iY, plint iZ) const
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
        return field[iX][iY][iZ];
    }
    T &operator[](plint ind)
    {
        PLB_PRECONDITION(ind >= 0 && ind < this->getNx() * this->getNy() * this->getNz());
        return rawData[ind];
    }
    T const &operator[](plint ind) const
    {
        PLB_PRECONDITION(ind >= 0 && ind < this->getNx() * this->getNy() * this->getNz());
        return rawData[ind];
    }

private:
    void allocateMemory();
    void releaseMemory();
    plint allocatedMemory() const;

private:
    bool ownsMemory;
    T *rawData;
    T ***field;

public:
    template <typename U>
    friend class NTensorField3D;
};
template <typename T, int nDim>
class TensorField3D;

template <typename T, int nDim>
class TensorFieldDataTransfer3D : public BlockDataTransfer3D {
public:
    TensorFieldDataTransfer3D();
    virtual void setBlock(AtomicBlock3D &block);
    virtual void setConstBlock(AtomicBlock3D const &block);
    virtual TensorFieldDataTransfer3D<T, nDim> *clone() const;
    virtual plint staticCellSize() const;
    /// Send data from the block into a byte-stream.
    virtual void send(Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the block.
    virtual void receive(Box3D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset)
    {
        receive(domain, buffer, kind);
    }
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    {
        receive(domain, buffer, kind);
    }
    /// Attribute data between two blocks.
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind);
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind, Dot3D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
    }

private:
    TensorField3D<T, nDim> *field;
    TensorField3D<T, nDim> const *constField;
};

template <typename T, int nDim>
class TensorField3D : public TensorFieldBase3D<T, nDim>, public AtomicBlock3D {
public:
    TensorField3D(plint nx_, plint ny_, plint nz_);
    TensorField3D(plint nx_, plint ny_, plint nz_, Array<T, nDim> const &iniVal);
    ~TensorField3D();
    TensorField3D(TensorField3D<T, nDim> const &rhs);
    TensorField3D(NTensorField3D<T> &rhs);
    TensorField3D<T, nDim> &operator=(TensorField3D<T, nDim> const &rhs);
    void swap(TensorField3D<T, nDim> &rhs);

public:
    virtual void reset();
    virtual Array<T, nDim> &get(plint iX, plint iY, plint iZ)
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
        return field[iX][iY][iZ];
    }
    virtual Array<T, nDim> const &get(plint iX, plint iY, plint iZ) const
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
        return field[iX][iY][iZ];
    }
    Array<T, nDim> &operator[](plint ind)
    {
        PLB_PRECONDITION(ind >= 0 && ind < this->getNx() * this->getNy() * this->getNz());
        return rawData[ind];
    }
    Array<T, nDim> const &operator[](plint ind) const
    {
        PLB_PRECONDITION(ind >= 0 && ind < this->getNx() * this->getNy() * this->getNz());
        return rawData[ind];
    }

private:
    void allocateMemory();
    void releaseMemory();
    plint allocatedMemory() const;

private:
    bool ownsMemory;
    Array<T, nDim> *rawData;
    Array<T, nDim> ***field;

public:
    template <typename U>
    friend class NTensorField3D;
};

template <typename T>
class NTensorField3D;

template <typename T>
class NTensorFieldDataTransfer3D : public BlockDataTransfer3D {
public:
    NTensorFieldDataTransfer3D();
    virtual void setBlock(AtomicBlock3D &block);
    virtual void setConstBlock(AtomicBlock3D const &block);
    virtual NTensorFieldDataTransfer3D<T> *clone() const;
    virtual plint staticCellSize() const;
    /// Send data from the block into a byte-stream.
    virtual void send(Box3D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the block.
    virtual void receive(Box3D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot3D absoluteOffset)
    {
        receive(domain, buffer, kind);
    }
    virtual void receive(
        Box3D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    {
        receive(domain, buffer, kind);
    }
    /// Attribute data between two blocks.
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind);
    virtual void attribute(
        Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ, AtomicBlock3D const &from,
        modif::ModifT kind, Dot3D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
    }

private:
    NTensorField3D<T> *field;
    NTensorField3D<T> const *constField;
};

template <typename T>
class NTensorField3D : public NTensorFieldBase3D<T>, public AtomicBlock3D {
public:
    NTensorField3D(plint nx_, plint ny_, plint nz_, plint ndim_);
    NTensorField3D(plint nx_, plint ny_, plint nz_, plint ndim_, T const *iniVal);
    ~NTensorField3D();
    NTensorField3D(NTensorField3D<T> const &rhs);
    NTensorField3D(ScalarField3D<T> &rhs);
    template <int nDim>
    NTensorField3D(TensorField3D<T, nDim> &rhs);
    NTensorField3D<T> &operator=(NTensorField3D<T> const &rhs);
    void swap(NTensorField3D<T> &rhs);

public:
    virtual void reset();
    virtual T *get(plint iX, plint iY, plint iZ)
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
        return field[iX][iY][iZ];
    }
    virtual T const *get(plint iX, plint iY, plint iZ) const
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        PLB_PRECONDITION(iZ >= 0 && iZ < this->getNz());
        return field[iX][iY][iZ];
    }
    T &operator[](plint ind)
    {
        PLB_PRECONDITION(
            ind >= 0 && ind < this->getNx() * this->getNy() * this->getNz() * this->getNdim());
        return rawData[ind];
    }
    T const &operator[](plint ind) const
    {
        PLB_PRECONDITION(
            ind >= 0 && ind < this->getNx() * this->getNy() * this->getNz() * this->getNdim());
        return rawData[ind];
    }

private:
    void allocateMemory();
    void releaseMemory();
    plint allocatedMemory() const;

private:
    bool ownsMemory;
    T *rawData;
    T ****field;

public:
    template <typename U>
    friend class ScalarField3D;
    template <typename U, int nDim>
    friend class TensorField3D;
};

}  // namespace plb

#endif  // DATA_FIELD_3D_H
