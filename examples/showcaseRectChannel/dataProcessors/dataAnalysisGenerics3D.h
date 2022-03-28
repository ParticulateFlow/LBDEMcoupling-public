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
 * Generic algorithms for data analysis.
 */
#ifndef DATA_ANALYSIS_GENERICS_3D_H
#define DATA_ANALYSIS_GENERICS_3D_H

#include "core/cell.h"

namespace plb {

/* ******** CountLatticeElementsFunctional3D **************************** */

template <typename T, template <typename U> class Descriptor, class BoolMask>
CountLatticeElementsFunctional3D<T, Descriptor, BoolMask>::CountLatticeElementsFunctional3D(
    BoolMask boolMask_) :
    countId(this->getStatistics().subscribeIntSum()), boolMask(boolMask_)
{ }

template <typename T, template <typename U> class Descriptor, class BoolMask>
void CountLatticeElementsFunctional3D<T, Descriptor, BoolMask>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                if (boolMask(cell)) {
                    statistics.gatherIntSum(countId, 1);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor, class BoolMask>
CountLatticeElementsFunctional3D<T, Descriptor, BoolMask>
    *CountLatticeElementsFunctional3D<T, Descriptor, BoolMask>::clone() const
{
    return new CountLatticeElementsFunctional3D<T, Descriptor, BoolMask>(*this);
}

template <typename T, template <typename U> class Descriptor, class BoolMask>
plint CountLatticeElementsFunctional3D<T, Descriptor, BoolMask>::getCount() const
{
    return this->getStatistics().getIntSum(countId);
}

/* ******** CountScalarElementsFunctional3D **************************** */

template <typename T, class BoolMask>
CountScalarElementsFunctional3D<T, BoolMask>::CountScalarElementsFunctional3D(BoolMask boolMask_) :
    countId(this->getStatistics().subscribeIntSum()), boolMask(boolMask_)
{ }

template <typename T, class BoolMask>
void CountScalarElementsFunctional3D<T, BoolMask>::process(Box3D domain, ScalarField3D<T> &field)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T value = field.get(iX, iY, iZ);
                if (boolMask(value)) {
                    statistics.gatherIntSum(countId, 1);
                }
            }
        }
    }
}

template <typename T, class BoolMask>
CountScalarElementsFunctional3D<T, BoolMask> *CountScalarElementsFunctional3D<T, BoolMask>::clone()
    const
{
    return new CountScalarElementsFunctional3D<T, BoolMask>(*this);
}

template <typename T, class BoolMask>
plint CountScalarElementsFunctional3D<T, BoolMask>::getCount() const
{
    return this->getStatistics().getIntSum(countId);
}

/* ******** CountTensorElementsFunctional3D **************************** */

template <typename T, int nDim, class BoolMask>
CountTensorElementsFunctional3D<T, nDim, BoolMask>::CountTensorElementsFunctional3D(
    BoolMask boolMask_) :
    countId(this->getStatistics().subscribeIntSum()), boolMask(boolMask_)
{ }

template <typename T, int nDim, class BoolMask>
void CountTensorElementsFunctional3D<T, nDim, BoolMask>::process(
    Box3D domain, TensorField3D<T, nDim> &field)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, nDim> const &value = field.get(iX, iY, iZ);
                if (boolMask(value)) {
                    statistics.gatherIntSum(countId, 1);
                }
            }
        }
    }
}

template <typename T, int nDim, class BoolMask>
CountTensorElementsFunctional3D<T, nDim, BoolMask>
    *CountTensorElementsFunctional3D<T, nDim, BoolMask>::clone() const
{
    return new CountTensorElementsFunctional3D<T, nDim, BoolMask>(*this);
}

template <typename T, int nDim, class BoolMask>
plint CountTensorElementsFunctional3D<T, nDim, BoolMask>::getCount() const
{
    return this->getStatistics().getIntSum(countId);
}

/* ******** ApplyScalarFunctional3D ************************************* */

template <typename T, class Function>
ApplyScalarFunctional3D<T, Function>::ApplyScalarFunctional3D(Function f_) : f(f_)
{ }

template <typename T, class Function>
void ApplyScalarFunctional3D<T, Function>::process(Box3D domain, ScalarField3D<T> &field)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                field.get(iX, iY, iZ) = f(field.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T, class Function>
ApplyScalarFunctional3D<T, Function> *ApplyScalarFunctional3D<T, Function>::clone() const
{
    return new ApplyScalarFunctional3D<T, Function>(*this);
}

template <typename T, class Function>
void ApplyScalarFunctional3D<T, Function>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, class Function>
BlockDomain::DomainT ApplyScalarFunctional3D<T, Function>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** EvaluateScalarFunctional3D ************************************* */

template <typename T, class Function>
EvaluateScalarFunctional3D<T, Function>::EvaluateScalarFunctional3D(Function f_) : f(f_)
{ }

template <typename T, class Function>
void EvaluateScalarFunctional3D<T, Function>::process(
    Box3D domain, ScalarField3D<T> &field, ScalarField3D<T> &result)
{
    Dot3D offset = computeRelativeDisplacement(field, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                result.get(iX + offset.x, iY + offset.y, iZ + offset.z) = f(field.get(iX, iY, iZ));
            }
        }
    }
}

template <typename T, class Function>
EvaluateScalarFunctional3D<T, Function> *EvaluateScalarFunctional3D<T, Function>::clone() const
{
    return new EvaluateScalarFunctional3D<T, Function>(*this);
}

template <typename T, class Function>
void EvaluateScalarFunctional3D<T, Function>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, class Function>
BlockDomain::DomainT EvaluateScalarFunctional3D<T, Function>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // DATA_ANALYSIS_GENERICS_3D_H
