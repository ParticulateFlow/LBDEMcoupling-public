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
 * Serial implementation of scalar, vector and tensor fields for 2D data analysis.
 * -- header file
 */

#ifndef DATA_FIELD_2D_H
#define DATA_FIELD_2D_H

#include "atomicBlock/atomicBlock2D.h"
#include "core/dataFieldBase2D.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"

namespace plb {

template <typename T>
class ScalarField2D;

template <typename T>
class ScalarFieldDataTransfer2D : public BlockDataTransfer2D {
public:
    ScalarFieldDataTransfer2D(ScalarField2D<T> &field_);
    virtual plint staticCellSize() const;
    /// Send data from the block into a byte-stream.
    virtual void send(Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the block.
    virtual void receive(Box2D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot2D offset)
    {
        receive(domain, buffer, kind);
    }
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    {
        receive(domain, buffer, kind);
    }
    /// Attribute data between two blocks.
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind);
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind,
        Dot2D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, from, kind);
    }

private:
    ScalarField2D<T> &field;
};

template <typename T>
class ScalarField2D : public ScalarFieldBase2D<T>, public AtomicBlock2D {
public:
    ScalarField2D(plint nx_, plint ny_, T iniVal = T());
    ~ScalarField2D();
    ScalarField2D(ScalarField2D<T> const &rhs);
    ScalarField2D<T> &operator=(ScalarField2D<T> const &rhs);
    void swap(ScalarField2D<T> &rhs);

public:
    virtual void reset();
    virtual pluint getSize() const
    {
        return (pluint)this->getNx() * (pluint)this->getNy();
    }
    virtual T &get(plint iX, plint iY)
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        return field[iX][iY];
    }
    virtual T const &get(plint iX, plint iY) const
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        return field[iX][iY];
    }
    T &operator[](plint ind)
    {
        PLB_PRECONDITION(ind >= 0 && ind < this->getNx() * this->getNy());
        return rawData[ind];
    }
    T const &operator[](plint ind) const
    {
        PLB_PRECONDITION(ind >= 0 && ind < this->getNx() * this->getNy());
        return rawData[ind];
    }
    /// Get access to data transfer between blocks
    virtual BlockDataTransfer2D &getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual BlockDataTransfer2D const &getDataTransfer() const;

private:
    void allocateMemory();
    void releaseMemory();

private:
    T *rawData;
    T **field;
    ScalarFieldDataTransfer2D<T> dataTransfer;
};

template <typename T, int nDim>
class TensorField2D;

template <typename T, int nDim>
class TensorFieldDataTransfer2D : public BlockDataTransfer2D {
public:
    TensorFieldDataTransfer2D(TensorField2D<T, nDim> &field_);
    virtual plint staticCellSize() const;
    /// Send data from the block into a byte-stream.
    virtual void send(Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the block.
    virtual void receive(Box2D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot2D offset)
    {
        receive(domain, buffer, kind);
    }
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    {
        return receive(domain, buffer, kind);
    }
    /// Attribute data between two blocks.
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind);
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind,
        Dot2D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, from, kind);
    }

private:
    TensorField2D<T, nDim> &field;
};

template <typename T, int nDim>
class TensorField2D : public TensorFieldBase2D<T, nDim>, public AtomicBlock2D {
public:
    TensorField2D(plint nx_, plint ny_);
    TensorField2D(plint nx_, plint ny_, Array<T, nDim> const &iniVal);
    ~TensorField2D();
    TensorField2D(TensorField2D<T, nDim> const &rhs);
    TensorField2D<T, nDim> &operator=(TensorField2D<T, nDim> const &rhs);
    void swap(TensorField2D<T, nDim> &rhs);

public:
    virtual void reset();
    virtual Array<T, nDim> &get(plint iX, plint iY)
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        return field[iX][iY];
    }
    virtual Array<T, nDim> const &get(plint iX, plint iY) const
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        return field[iX][iY];
    }
    Array<T, nDim> &operator[](plint ind)
    {
        PLB_PRECONDITION(ind >= 0 && ind < this->getNx() * this->getNy());
        return rawData[ind];
    }
    Array<T, nDim> const &operator[](plint ind) const
    {
        PLB_PRECONDITION(ind >= 0 && ind < this->getNx() * this->getNy());
        return rawData[ind];
    }
    /// Get access to data transfer between blocks
    virtual BlockDataTransfer2D &getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual BlockDataTransfer2D const &getDataTransfer() const;

private:
    void allocateMemory();
    void releaseMemory();

private:
    Array<T, nDim> *rawData;
    Array<T, nDim> **field;
    TensorFieldDataTransfer2D<T, nDim> dataTransfer;
};

template <typename T>
class NTensorField2D;

template <typename T>
class NTensorFieldDataTransfer2D : public BlockDataTransfer2D {
public:
    NTensorFieldDataTransfer2D(NTensorField2D<T> &field_);
    virtual plint staticCellSize() const;
    /// Send data from the block into a byte-stream.
    virtual void send(Box2D domain, std::vector<char> &buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the block.
    virtual void receive(Box2D domain, std::vector<char> const &buffer, modif::ModifT kind);
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind, Dot2D offset)
    {
        receive(domain, buffer, kind);
    }
    virtual void receive(
        Box2D domain, std::vector<char> const &buffer, modif::ModifT kind,
        std::map<int, std::string> const &foreignIds)
    {
        receive(domain, buffer, kind);
    }
    /// Attribute data between two blocks.
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind);
    virtual void attribute(
        Box2D toDomain, plint deltaX, plint deltaY, AtomicBlock2D const &from, modif::ModifT kind,
        Dot2D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, from, kind);
    }

private:
    NTensorField2D<T> &field;
};

template <typename T>
class NTensorField2D : public NTensorFieldBase2D<T>, public AtomicBlock2D {
public:
    NTensorField2D(plint nx_, plint ny_, plint ndim_);
    NTensorField2D(plint nx_, plint ny_, plint ndim_, T const *iniVal);
    ~NTensorField2D();
    NTensorField2D(NTensorField2D<T> const &rhs);
    NTensorField2D<T> &operator=(NTensorField2D<T> const &rhs);
    void swap(NTensorField2D<T> &rhs);

public:
    virtual void reset();
    virtual T *get(plint iX, plint iY)
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        return field[iX][iY];
    }
    virtual T const *get(plint iX, plint iY) const
    {
        PLB_PRECONDITION(iX >= 0 && iX < this->getNx());
        PLB_PRECONDITION(iY >= 0 && iY < this->getNy());
        return field[iX][iY];
    }
    T &operator[](plint ind)
    {
        PLB_PRECONDITION(ind >= 0 && ind < this->getNx() * this->getNy() * this->getNdim());
        return rawData[ind];
    }
    T const &operator[](plint ind) const
    {
        PLB_PRECONDITION(ind >= 0 && ind < this->getNx() * this->getNy() * this->getNdim());
        return rawData[ind];
    }
    /// Get access to data transfer between blocks
    virtual BlockDataTransfer2D &getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual BlockDataTransfer2D const &getDataTransfer() const;

private:
    void allocateMemory();
    void releaseMemory();

private:
    T *rawData;
    T ***field;
    NTensorFieldDataTransfer2D<T> dataTransfer;
};

}  // namespace plb

#endif  // DATA_FIELD_2D_H
