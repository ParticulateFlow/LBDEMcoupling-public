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
 * Scalar, vector and tensor fields for 2D data fields -- header file.
 */

#ifndef MULTI_DATA_FIELD_2D_H
#define MULTI_DATA_FIELD_2D_H

#include <vector>

#include "atomicBlock/dataField2D.h"
#include "core/dataFieldBase2D.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "multiBlock/multiBlock2D.h"

namespace plb {

template <typename T>
class MultiScalarField2D;

template <typename T>
struct MultiScalarAccess2D {
    virtual ~MultiScalarAccess2D() { }
    virtual T &getDistributedScalar(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, ScalarField2D<T> *> &fields) = 0;
    virtual T const &getDistributedScalar(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, ScalarField2D<T> *> const &fields) const = 0;
    virtual MultiScalarAccess2D<T> *clone() const = 0;
};

template <typename T>
class MultiScalarField2D : public ScalarFieldBase2D<T>, public MultiBlock2D {
public:
    typedef std::map<plint, ScalarField2D<T> *> BlockMap;

public:
    MultiScalarField2D(
        MultiBlockManagement2D const &multiBlockManagement_,
        BlockCommunicator2D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiScalarAccess2D<T> *multiScalarAccess_, T iniVal = T());
    MultiScalarField2D(plint nx, plint ny, T iniVal = T());
    ~MultiScalarField2D();
    MultiScalarField2D(MultiScalarField2D<T> const &rhs);
    MultiScalarField2D(MultiBlock2D const &rhs);
    /// Extract sub-domain from rhs and construct a multi-scalar-field with the same
    ///  data distribution and policy-classes; but the data itself and the data-processors
    ///  are not copied. MultiScalarAccess takes default value.
    MultiScalarField2D(MultiBlock2D const &rhs, Box2D subDomain, bool crop = true);
    MultiScalarField2D<T> &operator=(MultiScalarField2D<T> const &rhs);
    virtual MultiScalarField2D<T> *clone() const;
    virtual MultiScalarField2D<T> *clone(
        MultiBlockManagement2D const &newMultiBlockManagement) const;
    void swap(MultiScalarField2D<T> &rhs);

public:
    virtual void reset();
    virtual T &get(plint iX, plint iY);
    virtual T const &get(plint iX, plint iY) const;

public:
    virtual ScalarField2D<T> &getComponent(plint iBlock);
    virtual ScalarField2D<T> const &getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);
    virtual std::string getBlockName() const;
    virtual std::vector<std::string> getTypeInfo() const;
    static std::string blockName();
    static std::string basicType();

private:
    void allocateFields(T iniVal = T());
    void deAllocateFields();

private:
    BlockMap fields;
    MultiScalarAccess2D<T> *multiScalarAccess;

public:
    static const int staticId;
};

template <typename T>
MultiScalarField2D<T> &findMultiScalarField2D(id_t id);

template <typename T, int nDim>
class MultiTensorField2D;

template <typename T, int nDim>
struct MultiTensorAccess2D {
    virtual ~MultiTensorAccess2D() { }
    virtual Array<T, nDim> &getDistributedTensor(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, TensorField2D<T, nDim> *> &fields) = 0;
    virtual Array<T, nDim> const &getDistributedTensor(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, TensorField2D<T, nDim> *> const &fields) const = 0;
    virtual MultiTensorAccess2D<T, nDim> *clone() const = 0;
};

template <typename T, int nDim>
class MultiTensorField2D : public TensorFieldBase2D<T, nDim>, public MultiBlock2D {
public:
    typedef std::map<plint, TensorField2D<T, nDim> *> BlockMap;

public:
    MultiTensorField2D(
        MultiBlockManagement2D const &multiBlockManagement_,
        BlockCommunicator2D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiTensorAccess2D<T, nDim> *multiTensorAccess_);
    MultiTensorField2D(
        MultiBlockManagement2D const &multiBlockManagement_,
        BlockCommunicator2D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiTensorAccess2D<T, nDim> *multiTensorAccess_, Array<T, nDim> const &iniVal);
    MultiTensorField2D(plint nx, plint ny);
    MultiTensorField2D(plint nx, plint ny, Array<T, nDim> const &iniVal);
    ~MultiTensorField2D();
    MultiTensorField2D(MultiTensorField2D<T, nDim> const &rhs);
    MultiTensorField2D(MultiBlock2D const &rhs);
    /// Extract sub-domain from rhs and construct a multi-tensor-field with the same
    ///  data distribution and policy-classes; but the data itself and the data-processors
    ///  are not copied. MultiTensorAccess takes default value.
    MultiTensorField2D(MultiBlock2D const &rhs, Box2D subDomain, bool crop = true);
    MultiTensorField2D<T, nDim> &operator=(MultiTensorField2D<T, nDim> const &rhs);
    virtual MultiTensorField2D<T, nDim> *clone() const;
    virtual MultiTensorField2D<T, nDim> *clone(
        MultiBlockManagement2D const &newMultiBlockManagement) const;
    void swap(MultiTensorField2D<T, nDim> &rhs);

public:
    virtual void reset();
    virtual Array<T, nDim> &get(plint iX, plint iY);
    virtual Array<T, nDim> const &get(plint iX, plint iY) const;

public:
    virtual TensorField2D<T, nDim> &getComponent(plint iBlock);
    virtual TensorField2D<T, nDim> const &getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);
    virtual std::string getBlockName() const;
    virtual std::vector<std::string> getTypeInfo() const;
    static std::string blockName();
    static std::string basicType();

private:
    void allocateFields();
    void allocateFields(Array<T, nDim> const &iniVal);
    void deAllocateFields();

private:
    BlockMap fields;
    MultiTensorAccess2D<T, nDim> *multiTensorAccess;

public:
    static const int staticId;
};

template <typename T, int nDim>
MultiTensorField2D<T, nDim> &findMultiTensorField2D(id_t id);

template <typename T>
class MultiNTensorField2D;

template <typename T>
struct MultiNTensorAccess2D {
    virtual ~MultiNTensorAccess2D() { }
    virtual T *getDistributedNTensor(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, NTensorField2D<T> *> &fields) = 0;
    virtual T const *getDistributedNTensor(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, NTensorField2D<T> *> const &fields) const = 0;
    virtual MultiNTensorAccess2D<T> *clone() const = 0;
};

template <typename T>
class MultiNTensorField2D : public NTensorFieldBase2D<T>, public MultiBlock2D {
public:
    typedef std::map<plint, NTensorField2D<T> *> BlockMap;

public:
    MultiNTensorField2D(
        plint ndim, MultiBlockManagement2D const &multiBlockManagement_,
        BlockCommunicator2D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiNTensorAccess2D<T> *multiNTensorAccess_);
    MultiNTensorField2D(
        plint ndim, T const *iniVal, MultiBlockManagement2D const &multiBlockManagement_,
        BlockCommunicator2D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiNTensorAccess2D<T> *multiNTensorAccess_);
    MultiNTensorField2D(plint nx, plint ny, plint ndim);
    MultiNTensorField2D(plint nx, plint ny, plint ndim, T const *iniVal);
    ~MultiNTensorField2D();
    MultiNTensorField2D(MultiNTensorField2D<T> const &rhs);
    MultiNTensorField2D(plint ndim, MultiBlock2D const &rhs);
    /// Extract sub-domain from rhs and construct a multi-tensor-field with the same
    ///  data distribution and policy-classes; but the data itself and the data-processors
    ///  are not copied. MultiNTensorAccess takes default value.
    MultiNTensorField2D(plint ndim, MultiBlock2D const &rhs, Box2D subDomain, bool crop = true);
    MultiNTensorField2D<T> &operator=(MultiNTensorField2D<T> const &rhs);
    virtual MultiNTensorField2D<T> *clone() const;
    virtual MultiNTensorField2D<T> *clone(
        MultiBlockManagement2D const &newMultiBlockManagement) const;
    void swap(MultiNTensorField2D<T> &rhs);

public:
    virtual void reset();
    virtual T *get(plint iX, plint iY);
    virtual T const *get(plint iX, plint iY) const;

public:
    virtual NTensorField2D<T> &getComponent(plint iBlock);
    virtual NTensorField2D<T> const &getComponent(plint iBlock) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock2D const &fromBlock, Box2D const &fromDomain, Box2D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);
    virtual std::string getBlockName() const;
    virtual std::vector<std::string> getTypeInfo() const;
    static std::string blockName();
    static std::string basicType();

private:
    void allocateFields();
    void allocateFields(T const *iniVal);
    void deAllocateFields();

private:
    BlockMap fields;
    MultiNTensorAccess2D<T> *multiNTensorAccess;

public:
    static const int staticId;
};

template <typename T>
MultiNTensorField2D<T> &findMultiNTensorField2D(id_t id);

}  // namespace plb

#endif  // MULTI_DATA_FIELD_2D_H
