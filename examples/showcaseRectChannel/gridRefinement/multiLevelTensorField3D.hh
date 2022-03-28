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

#ifndef MULTI_LEVEL_TENSOR_FIELD_3D_HH
#define MULTI_LEVEL_TENSOR_FIELD_3D_HH

#include <map>
#include <vector>

#include "core/globalDefs.h"
#include "gridRefinement/boxLogic3D.hh"
#include "gridRefinement/multiLevelTensorField3D.h"
#include "gridRefinement/octreeGridStructure.h"

namespace plb {

// ======================================================================= //
// ====================MultiLevelTensorField3D============================ //
// ======================================================================= //

template <typename T, int nDim>
MultiLevelTensorField3D<T, nDim>::MultiLevelTensorField3D(const OctreeGridStructure &ogs_) :
    MultiLevel3D(), ogs(ogs_)
{
    for (plint iLevel = 0; iLevel < getNumLevels(); ++iLevel) {
        gridLevels.push_back(new MultiTensorField3D<T, nDim>(
            ogs.getMultiBlockManagement(iLevel, 1),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));
    }
}

template <typename T, int nDim>
MultiLevelTensorField3D<T, nDim>::MultiLevelTensorField3D(
    const OctreeGridStructure &ogs_, Box3D domain, plint domainLevel) :
    MultiLevel3D(), ogs(ogs_)
{
    for (plint iLevel = 0; iLevel < getNumLevels(); ++iLevel) {
        Box3D currentDomain = (iLevel - domainLevel >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iLevel - domainLevel))
                                  : domain.divide(util::intTwoToThePower(domainLevel - iLevel));

        gridLevels.push_back(new MultiTensorField3D<T, nDim>(
            ogs.getMultiBlockManagement(iLevel, currentDomain, 1),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));
    }
}

template <typename T, int nDim>
MultiLevelTensorField3D<T, nDim>::MultiLevelTensorField3D(
    const OctreeGridStructure &ogs_, const Array<T, nDim> &iniVal) :
    MultiLevel3D(), ogs(ogs_)
{
    for (plint iLevel = 0; iLevel < getNumLevels(); ++iLevel) {
        gridLevels.push_back(new MultiTensorField3D<T, nDim>(
            ogs.getMultiBlockManagement(iLevel, 1),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>(), iniVal));
    }
}

template <typename T, int nDim>
MultiLevelTensorField3D<T, nDim>::MultiLevelTensorField3D(
    const OctreeGridStructure &ogs_, Box3D domain, plint domainLevel,
    const Array<T, nDim> &iniVal) :
    MultiLevel3D(), ogs(ogs_)
{
    for (plint iLevel = 0; iLevel < getNumLevels(); ++iLevel) {
        Box3D currentDomain = (iLevel - domainLevel >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iLevel - domainLevel))
                                  : domain.divide(util::intTwoToThePower(domainLevel - iLevel));

        gridLevels.push_back(new MultiTensorField3D<T, nDim>(
            ogs.getMultiBlockManagement(iLevel, currentDomain, 1),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>(), iniVal));
    }
}

template <typename T, int nDim>
MultiLevelTensorField3D<T, nDim>::MultiLevelTensorField3D(
    MultiLevelTensorField3D<T, nDim> const &rhs) :
    MultiLevel3D(), ogs(rhs.ogs), gridLevels(rhs.gridLevels)
{ }

template <typename T, int nDim>
MultiLevelTensorField3D<T, nDim> &MultiLevelTensorField3D<T, nDim>::operator=(
    MultiLevelTensorField3D<T, nDim> const &rhs)
{
    MultiLevelTensorField3D<T, nDim>(rhs).swap(*this);
    return *this;
}

template <typename T, int nDim>
void MultiLevelTensorField3D<T, nDim>::swap(MultiLevelTensorField3D<T, nDim> &rhs)
{
    std::swap(ogs, rhs.ogs);
    std::swap(gridLevels, rhs.gridLevels);
}

template <typename T, int nDim>
MultiLevelTensorField3D<T, nDim>::~MultiLevelTensorField3D()
{
    for (pluint iA = 0; iA < gridLevels.size(); ++iA) {
        delete gridLevels[iA];
    }
}

// ======================================================================= //
// ====================MultiLevelTensorFieldForOutput3D=================== //
// ======================================================================= //

template <typename T, int nDim>
MultiLevelTensorFieldForOutput3D<T, nDim>::MultiLevelTensorFieldForOutput3D(
    const OctreeGridStructure &ogs_, bool crop) :
    MultiLevel3D(), ogs(ogs_)
{
    for (plint iLevel = 0; iLevel < getNumLevels(); ++iLevel) {
        gridLevels.push_back(new MultiTensorField3D<T, nDim>(
            ogs.getMultiBlockManagementForOutput(iLevel, crop, 1),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));
    }
}

template <typename T, int nDim>
MultiLevelTensorFieldForOutput3D<T, nDim>::MultiLevelTensorFieldForOutput3D(
    const OctreeGridStructure &ogs_, Box3D domain, plint domainLevel, bool crop) :
    MultiLevel3D(), ogs(ogs_)
{
    for (plint iLevel = 0; iLevel < getNumLevels(); ++iLevel) {
        Box3D currentDomain = (iLevel - domainLevel >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iLevel - domainLevel))
                                  : domain.divide(util::intTwoToThePower(domainLevel - iLevel));

        gridLevels.push_back(new MultiTensorField3D<T, nDim>(
            ogs.getMultiBlockManagementForOutput(iLevel, currentDomain, crop, 1),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiTensorAccess<T, nDim>()));
    }
}

template <typename T, int nDim>
MultiLevelTensorFieldForOutput3D<T, nDim>::MultiLevelTensorFieldForOutput3D(
    MultiLevelTensorFieldForOutput3D<T, nDim> const &rhs) :
    MultiLevel3D(), ogs(rhs.ogs), gridLevels(rhs.gridLevels)
{ }

template <typename T, int nDim>
MultiLevelTensorFieldForOutput3D<T, nDim> &MultiLevelTensorFieldForOutput3D<T, nDim>::operator=(
    MultiLevelTensorFieldForOutput3D<T, nDim> const &rhs)
{
    MultiLevelTensorFieldForOutput3D<T, nDim>(rhs).swap(*this);
    return *this;
}

template <typename T, int nDim>
void MultiLevelTensorFieldForOutput3D<T, nDim>::swap(MultiLevelTensorFieldForOutput3D<T, nDim> &rhs)
{
    std::swap(ogs, rhs.ogs);
    std::swap(gridLevels, rhs.gridLevels);
}

template <typename T, int nDim>
MultiLevelTensorFieldForOutput3D<T, nDim>::~MultiLevelTensorFieldForOutput3D()
{
    for (pluint iA = 0; iA < gridLevels.size(); ++iA) {
        delete gridLevels[iA];
    }
}

}  // namespace plb

#endif  // MULTI_LEVEL_TENSOR_FIELD_3D_HH
