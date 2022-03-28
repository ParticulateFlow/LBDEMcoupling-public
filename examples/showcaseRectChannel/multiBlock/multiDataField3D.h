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
 * Scalar, vector and tensor fields for 3D data fields -- header file.
 */

#ifndef MULTI_DATA_FIELD_3D_H
#define MULTI_DATA_FIELD_3D_H

#include <vector>

#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataField3D.h"
#include "core/dataFieldBase2D.h"
#include "core/dataFieldBase3D.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "multiBlock/multiBlock3D.h"

namespace plb {

template <typename T>
class MultiScalarField3D;
template <typename T>
class MultiNTensorField3D;
template <typename T, int nDim>
class MultiTensorField3D;

template <typename T>
struct MultiScalarAccess3D {
    virtual ~MultiScalarAccess3D() { }
    virtual T &getDistributedScalar(
        plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
        std::map<plint, ScalarField3D<T> *> &fields) = 0;
    virtual T const &getDistributedScalar(
        plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
        std::map<plint, ScalarField3D<T> *> const &fields) const = 0;
    virtual MultiScalarAccess3D<T> *clone() const = 0;
};

template <typename T>
class MultiScalarField3D : public ScalarFieldBase3D<T>, public MultiBlock3D {
public:
    typedef std::map<plint, ScalarField3D<T> *> BlockMap;

public:
    MultiScalarField3D(
        MultiBlockManagement3D const &multiBlockManagement_,
        BlockCommunicator3D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiScalarAccess3D<T> *multiScalarAccess_, T iniVal = T());
    MultiScalarField3D(plint nx, plint ny, plint nz, T iniVal = T());
    ~MultiScalarField3D();
    MultiScalarField3D(MultiScalarField3D<T> const &rhs);
    MultiScalarField3D(MultiNTensorField3D<T> &rhs, bool shareMemory);
    MultiScalarField3D(MultiBlock3D const &rhs);
    /// Extract sub-domain from rhs and construct a multi-scalar-field with the same
    ///  data distribution and policy-classes; but the data itself and the data-processors
    ///  are not copied. MultiScalarAccess takes default value.
    MultiScalarField3D(MultiBlock3D const &rhs, Box3D subDomain, bool crop = true);
    MultiScalarField3D<T> &operator=(MultiScalarField3D<T> const &rhs);
    virtual MultiScalarField3D<T> *clone() const;
    virtual MultiScalarField3D<T> *clone(
        MultiBlockManagement3D const &newMultiBlockManagement) const;
    void swap(MultiScalarField3D<T> &rhs);

public:
    virtual void reset();
    virtual T &get(plint iX, plint iY, plint iZ);
    virtual T const &get(plint iX, plint iY, plint iZ) const;

public:
    virtual ScalarField3D<T> &getComponent(plint blockId);
    virtual ScalarField3D<T> const &getComponent(plint blockId) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);
    virtual std::string getBlockName() const;
    virtual std::vector<std::string> getTypeInfo() const;
    static std::string blockName();
    static std::string basicType();

public:
    MultiNTensorField3D<T> &nTensorView();

private:
    void allocateFields(T iniVal = T());
    void deAllocateFields();

private:
    BlockMap fields;
    MultiScalarAccess3D<T> *multiScalarAccess;
    MultiNTensorField3D<T> *nTensorViewBlock;

public:
    static const int staticId;
    template <typename U>
    friend class MultiNTensorField3D;
};

template <typename T>
MultiScalarField3D<T> &findMultiScalarField3D(id_t id);

template <typename T, int nDim>
struct MultiTensorAccess3D {
    virtual ~MultiTensorAccess3D() { }
    virtual Array<T, nDim> &getDistributedTensor(
        plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
        std::map<plint, TensorField3D<T, nDim> *> &fields) = 0;
    virtual Array<T, nDim> const &getDistributedTensor(
        plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
        std::map<plint, TensorField3D<T, nDim> *> const &fields) const = 0;
    virtual MultiTensorAccess3D<T, nDim> *clone() const = 0;
};

template <typename T, int nDim>
class MultiTensorField3D : public TensorFieldBase3D<T, nDim>, public MultiBlock3D {
public:
    typedef std::map<plint, TensorField3D<T, nDim> *> BlockMap;

public:
    MultiTensorField3D(
        MultiBlockManagement3D const &multiBlockManagement_,
        BlockCommunicator3D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiTensorAccess3D<T, nDim> *multiTensorAccess_);
    MultiTensorField3D(
        MultiBlockManagement3D const &multiBlockManagement_,
        BlockCommunicator3D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiTensorAccess3D<T, nDim> *multiTensorAccess_, Array<T, nDim> const &iniVal);
    MultiTensorField3D(plint nx, plint ny, plint nz);
    MultiTensorField3D(plint nx, plint ny, plint nz, Array<T, nDim> const &iniVal);
    ~MultiTensorField3D();
    MultiTensorField3D(MultiTensorField3D<T, nDim> const &rhs);
    MultiTensorField3D(MultiNTensorField3D<T> &rhs, bool shareMemory);
    MultiTensorField3D(MultiBlock3D const &rhs);
    /// Extract sub-domain from rhs and construct a multi-tensor-field with the same
    ///  data distribution and policy-classes; but the data itself and the data-processors
    ///  are not copied. MultiTensorAccess takes default value.
    MultiTensorField3D(MultiBlock3D const &rhs, Box3D subDomain, bool crop = true);
    MultiTensorField3D<T, nDim> &operator=(MultiTensorField3D<T, nDim> const &rhs);
    virtual MultiTensorField3D<T, nDim> *clone() const;
    virtual MultiTensorField3D<T, nDim> *clone(
        MultiBlockManagement3D const &newMultiBlockManagement) const;
    void swap(MultiTensorField3D<T, nDim> &rhs);

public:
    virtual void reset();
    virtual Array<T, nDim> &get(plint iX, plint iY, plint iZ);
    virtual Array<T, nDim> const &get(plint iX, plint iY, plint iZ) const;

public:
    virtual TensorField3D<T, nDim> &getComponent(plint blockId);
    virtual TensorField3D<T, nDim> const &getComponent(plint blockId) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);
    virtual std::string getBlockName() const;
    virtual std::vector<std::string> getTypeInfo() const;
    static std::string blockName();
    static std::string basicType();

public:
    MultiNTensorField3D<T> &nTensorView();

private:
    void allocateFields();
    void allocateFields(Array<T, nDim> const &iniVal);
    void deAllocateFields();

private:
    BlockMap fields;
    MultiTensorAccess3D<T, nDim> *multiTensorAccess;
    MultiNTensorField3D<T> *nTensorViewBlock;

public:
    static const int staticId;
    template <typename U>
    friend class MultiNTensorField3D;
};

template <typename T, int nDim>
MultiTensorField3D<T, nDim> &findMultiTensorField3D(id_t id);

template <typename T>
struct MultiNTensorAccess3D {
    virtual ~MultiNTensorAccess3D() { }
    virtual T *getDistributedNTensor(
        plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
        std::map<plint, NTensorField3D<T> *> &fields) = 0;
    virtual T const *getDistributedNTensor(
        plint iX, plint iY, plint iZ, MultiBlockManagement3D const &multiBlockManagement,
        std::map<plint, NTensorField3D<T> *> const &fields) const = 0;
    virtual MultiNTensorAccess3D<T> *clone() const = 0;
};

template <typename T>
class MultiNTensorField3D : public NTensorFieldBase3D<T>, public MultiBlock3D {
public:
    typedef std::map<plint, NTensorField3D<T> *> BlockMap;

public:
    MultiNTensorField3D(
        plint ndim, MultiBlockManagement3D const &multiBlockManagement_,
        BlockCommunicator3D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiNTensorAccess3D<T> *multiNTensorAccess_);
    MultiNTensorField3D(
        plint ndim, T const *iniVal, MultiBlockManagement3D const &multiBlockManagement_,
        BlockCommunicator3D *blockCommunicator_, CombinedStatistics *combinedStatistics_,
        MultiNTensorAccess3D<T> *multiNTensorAccess_);
    MultiNTensorField3D(plint nx, plint ny, plint nz, plint ndim);
    MultiNTensorField3D(plint nx, plint ny, plint nz, plint ndim, T const *iniVal);
    ~MultiNTensorField3D();
    MultiNTensorField3D(MultiNTensorField3D<T> const &rhs);
    MultiNTensorField3D(MultiScalarField3D<T> &rhs, bool shareMemory);
    template <int nDim>
    MultiNTensorField3D(MultiTensorField3D<T, nDim> &rhs, bool shareMemory);
    MultiNTensorField3D(plint ndim, MultiBlock3D const &rhs);
    /// Extract sub-domain from rhs and construct a multi-tensor-field with the same
    ///  data distribution and policy-classes; but the data itself and the data-processors
    ///  are not copied. MultiNTensorAccess takes default value.
    MultiNTensorField3D(plint ndim, MultiBlock3D const &rhs, Box3D subDomain, bool crop = true);
    MultiNTensorField3D<T> &operator=(MultiNTensorField3D<T> const &rhs);
    virtual MultiNTensorField3D<T> *clone() const;
    virtual MultiNTensorField3D<T> *clone(
        MultiBlockManagement3D const &newMultiBlockManagement) const;
    void swap(MultiNTensorField3D<T> &rhs);

public:
    virtual void reset();
    virtual T *get(plint iX, plint iY, plint iZ);
    virtual T const *get(plint iX, plint iY, plint iZ) const;

public:
    virtual NTensorField3D<T> &getComponent(plint blockId);
    virtual NTensorField3D<T> const &getComponent(plint blockId) const;
    virtual plint sizeOfCell() const;
    virtual plint getCellDim() const;
    virtual int getStaticId() const;
    virtual void copyReceive(
        MultiBlock3D const &fromBlock, Box3D const &fromDomain, Box3D const &toDomain,
        modif::ModifT whichData = modif::dataStructure);
    virtual std::string getBlockName() const;
    virtual std::vector<std::string> getTypeInfo() const;
    static std::string blockName();
    static std::string basicType();

public:
    MultiScalarField3D<T> &scalarView();
    template <int nDim>
    MultiTensorField3D<T, nDim> &tensorView();

private:
    void allocateFields();
    void allocateFields(T const *iniVal);
    void deAllocateFields();

private:
    BlockMap fields;
    MultiNTensorAccess3D<T> *multiNTensorAccess;
    MultiBlock3D *scalarOrTensorView;

public:
    template <typename U>
    friend class MultiScalarField3D;
    template <typename U, int nDim>
    friend class MultiTensorField3D;
    static const int staticId;
};

template <typename T>
MultiNTensorField3D<T> &findMultiNTensorField3D(id_t id);

}  // namespace plb

#endif  // MULTI_DATA_FIELD_3D_H
