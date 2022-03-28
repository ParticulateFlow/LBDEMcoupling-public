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
 * Coupling between grids of different refinement level -- header file.
 */

#ifndef MULTI_LEVEL_NTENSOR_FIELD_3D_HH
#define MULTI_LEVEL_NTENSOR_FIELD_3D_HH

#include <map>
#include <vector>

#include "core/globalDefs.h"
#include "gridRefinement/boxLogic3D.hh"
#include "gridRefinement/multiLevelNTensorField3D.h"
#include "gridRefinement/octreeGridStructure.h"

namespace plb {

// ======================================================================= //
// ====================MultiLevelNTensorField3D=========================== //
// ======================================================================= //

template <typename T>
MultiLevelNTensorField3D<T>::MultiLevelNTensorField3D(OctreeGridStructure &ogs_) :
    MultiLevel3D(), ogs(ogs_)
{
    for (plint iLevel = 0; iLevel < getNumLevels(); ++iLevel) {
        gridLevels.push_back(new MultiNTensorField3D<T>(
            ogs.getMultiBlockManagement(iLevel, 1),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>()));
    }
}

template <typename T>
MultiLevelNTensorField3D<T>::MultiLevelNTensorField3D(
    OctreeGridStructure &ogs_, Box3D domain, plint domainLevel) :
    MultiLevel3D(), ogs(ogs_)
{
    for (plint iLevel = 0; iLevel < getNumLevels(); ++iLevel) {
        Box3D currentDomain = domain.multiply(util::intTwoToThePower(iLevel - domainLevel));

        gridLevels.push_back(new MultiNTensorField3D<T>(
            ogs.getMultiBlockManagement(iLevel, currentDomain, 1),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>()));
    }
}

template <typename T>
MultiLevelNTensorField3D<T>::MultiLevelNTensorField3D(MultiLevelNTensorField3D<T> const &rhs) :
    MultiLevel3D(), ogs(rhs.ogs), gridLevels(rhs.gridLevels)
{ }

template <typename T>
MultiLevelNTensorField3D<T> &MultiLevelNTensorField3D<T>::operator=(
    MultiLevelNTensorField3D<T> const &rhs)
{
    MultiLevelNTensorField3D<T>(rhs).swap(*this);
    return *this;
}

template <typename T>
void MultiLevelNTensorField3D<T>::swap(MultiLevelNTensorField3D<T> &rhs)
{
    std::swap(ogs, rhs.ogs);
    std::swap(gridLevels, rhs.gridLevels);
}

template <typename T>
MultiLevelNTensorField3D<T>::~MultiLevelNTensorField3D()
{
    for (pluint iA = 0; iA < gridLevels.size(); ++iA) {
        delete gridLevels[iA];
    }
}

// ======================================================================= //
// ====================MultiLevelNTensorFieldForOutput3D=========================== //
// ======================================================================= //

template <typename T>
MultiLevelNTensorFieldForOutput3D<T>::MultiLevelNTensorFieldForOutput3D(
    OctreeGridStructure &ogs_, bool crop) :
    MultiLevel3D(), ogs(ogs_)
{
    for (plint iLevel = 0; iLevel < getNumLevels(); ++iLevel) {
        gridLevels.push_back(new MultiNTensorField3D<T>(
            ogs.getMultiBlockManagementForOutput(iLevel, crop, 1),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>()));
    }
}

template <typename T>
MultiLevelNTensorFieldForOutput3D<T>::MultiLevelNTensorFieldForOutput3D(
    OctreeGridStructure &ogs_, Box3D domain, plint domainLevel, bool crop) :
    MultiLevel3D(), ogs(ogs_)
{
    for (plint iLevel = 0; iLevel < getNumLevels(); ++iLevel) {
        Box3D currentDomain = domain.multiply(util::intTwoToThePower(iLevel - domainLevel));

        gridLevels.push_back(new MultiNTensorField3D<T>(
            ogs.getMultiBlockManagementForOutput(iLevel, currentDomain, crop, 1),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiNTensorAccess<T>()));
    }
}

template <typename T>
MultiLevelNTensorFieldForOutput3D<T>::MultiLevelNTensorFieldForOutput3D(
    MultiLevelNTensorFieldForOutput3D<T> const &rhs) :
    MultiLevel3D(), ogs(rhs.ogs), gridLevels(rhs.gridLevels)
{ }

template <typename T>
MultiLevelNTensorFieldForOutput3D<T> &MultiLevelNTensorFieldForOutput3D<T>::operator=(
    MultiLevelNTensorFieldForOutput3D<T> const &rhs)
{
    MultiLevelNTensorFieldForOutput3D<T>(rhs).swap(*this);
    return *this;
}

template <typename T>
void MultiLevelNTensorFieldForOutput3D<T>::swap(MultiLevelNTensorFieldForOutput3D<T> &rhs)
{
    std::swap(ogs, rhs.ogs);
    std::swap(gridLevels, rhs.gridLevels);
}

template <typename T>
MultiLevelNTensorFieldForOutput3D<T>::~MultiLevelNTensorFieldForOutput3D()
{
    for (pluint iA = 0; iA < gridLevels.size(); ++iA) {
        delete gridLevels[iA];
    }
}

}  // namespace plb

#endif  // MULTI_LEVEL_NTENSOR_FIELD_3D_HH
