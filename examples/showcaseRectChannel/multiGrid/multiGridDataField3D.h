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

#ifndef MULTI_GRID_DATA_FIELD_3D_H
#define MULTI_GRID_DATA_FIELD_3D_H

#include "atomicBlock/dataField3D.h"
#include "core/dataFieldBase3D.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "multiBlock/multiBlock3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiGrid/multiGrid3D.h"

namespace plb {

template <typename T>
class MultiGridScalarField3D : public ScalarFieldBase3D<T>, public MultiGrid3D {
public:
    MultiGridScalarField3D(
        MultiGridManagement3D management_, std::vector<BlockCommunicator3D *> communicators_,
        std::vector<CombinedStatistics *> combinedStatistics_, plint behaviorLevel_ = 0);

    MultiGridScalarField3D(MultiGridManagement3D management_, plint behaviorLevel_ = 0);

    MultiGridScalarField3D(MultiGridScalarField3D<T> const &rhs);
    MultiGridScalarField3D(MultiGrid3D const &rhs);
    MultiGridScalarField3D(MultiGrid3D const &rhs, Box3D subDomain, bool crop = true);

    ~MultiGridScalarField3D();

    /* ************ ScalarFieldBase3D ************* */
public:
    void reset();
    T &get(plint iX, plint iY, plint iZ);
    T const &get(plint iX, plint iY, plint iZ) const;

    /* ************ MultiGrid3D ************* */
public:
    MultiScalarField3D<T> &getComponent(plint level);
    MultiScalarField3D<T> const &getComponent(plint level) const;
    int getBlockId() const;
    std::unique_ptr<MultiScalarField3D<T> > convertToCoarsest(plint dimDx = 0, plint dimDt = 0);
    std::unique_ptr<MultiScalarField3D<T> > convertToFinest(plint dimDx = 0, plint dimDt = 0);
    std::unique_ptr<MultiScalarField3D<T> > convertToLevel(
        plint level, plint dimDx = 0, plint dimDt = 0);

private:
    /// Create the multiScalarFields for each level
    void allocateFields();
    void allocateFields(
        std::vector<BlockCommunicator3D *> communicators,
        std::vector<CombinedStatistics *> combinedStatistics);

private:
    std::vector<MultiScalarField3D<T> *> fields;
};

template <typename T, int nDim>
class MultiGridTensorField3D : public TensorFieldBase3D<T, nDim>, public MultiGrid3D {
public:
    MultiGridTensorField3D(
        MultiGridManagement3D management_, std::vector<BlockCommunicator3D *> communicators_,
        std::vector<CombinedStatistics *> combinedStatistics_, plint behaviorLevel_ = 0);

    MultiGridTensorField3D(MultiGridManagement3D management_, plint behaviorLevel_ = 0);

    MultiGridTensorField3D(MultiGridTensorField3D<T, nDim> const &rhs);
    MultiGridTensorField3D(MultiGrid3D const &rhs);
    MultiGridTensorField3D(MultiGrid3D const &rhs, Box3D subDomain, bool crop = true);

    ~MultiGridTensorField3D();

    /* ************ TensorFieldBase3D ************* */
public:
    void reset();
    virtual Array<T, nDim> &get(plint iX, plint iY, plint iZ);
    virtual Array<T, nDim> const &get(plint iX, plint iY, plint iZ) const;

    /* ************ MultiGrid3D ************* */
public:
    MultiTensorField3D<T, nDim> &getComponent(plint level);
    MultiTensorField3D<T, nDim> const &getComponent(plint level) const;
    int getBlockId() const;
    std::unique_ptr<MultiTensorField3D<T, nDim> > convertToCoarsest(
        plint dimDx = 0, plint dimDt = 0);
    std::unique_ptr<MultiTensorField3D<T, nDim> > convertToFinest(plint dimDx = 0, plint dimDt = 0);
    std::unique_ptr<MultiTensorField3D<T, nDim> > convertToLevel(
        plint level, plint dimDx = 0, plint dimDt = 0);

private:
    /// Create the multiTensorFields for each level
    void allocateFields();
    void allocateFields(
        std::vector<BlockCommunicator3D *> communicators,
        std::vector<CombinedStatistics *> combinedStatistics);

private:
    std::vector<MultiTensorField3D<T, nDim> *> fields;
};

}  // namespace plb

#endif  // MULTI_GRID_DATA_FIELD_3D_H
