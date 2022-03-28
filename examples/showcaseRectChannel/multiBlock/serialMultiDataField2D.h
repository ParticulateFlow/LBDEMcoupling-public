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
 * Serial access to elements of a scalar/tensor field -- header file.
 */
#ifndef SERIAL_MULTI_DATA_FIELD_2D_H
#define SERIAL_MULTI_DATA_FIELD_2D_H

#include "core/globalDefs.h"
#include "multiBlock/multiDataField2D.h"

namespace plb {

template <typename T>
class SerialScalarAccess2D : public MultiScalarAccess2D<T> {
public:
    SerialScalarAccess2D();
    virtual T &getDistributedScalar(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, ScalarField2D<T> *> &fields);
    virtual T const &getDistributedScalar(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, ScalarField2D<T> *> const &fields) const;
    virtual SerialScalarAccess2D<T> *clone() const;

private:
    mutable plint locatedBlock;
};

template <typename T, int nDim>
class SerialTensorAccess2D : public MultiTensorAccess2D<T, nDim> {
public:
    SerialTensorAccess2D();
    virtual Array<T, nDim> &getDistributedTensor(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, TensorField2D<T, nDim> *> &fields);
    virtual Array<T, nDim> const &getDistributedTensor(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, TensorField2D<T, nDim> *> const &fields) const;
    virtual SerialTensorAccess2D<T, nDim> *clone() const;

private:
    mutable plint locatedBlock;
};

template <typename T>
class SerialNTensorAccess2D : public MultiNTensorAccess2D<T> {
public:
    SerialNTensorAccess2D();
    virtual T *getDistributedNTensor(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, NTensorField2D<T> *> &fields);
    virtual T const *getDistributedNTensor(
        plint iX, plint iY, MultiBlockManagement2D const &multiBlockManagement,
        std::map<plint, NTensorField2D<T> *> const &fields) const;
    virtual SerialNTensorAccess2D<T> *clone() const;

private:
    mutable plint locatedBlock;
};

}  // namespace plb

#endif  // SERIAL_MULTI_DATA_FIELD_2D_H
