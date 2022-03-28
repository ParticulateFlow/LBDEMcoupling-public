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

#ifndef MULTI_LEVEL_TENSOR_FIELD_3D_H
#define MULTI_LEVEL_TENSOR_FIELD_3D_H

#include <map>
#include <vector>

#include "core/globalDefs.h"
#include "gridRefinement/boxLogic3D.h"
#include "gridRefinement/multiLevel3D.h"
#include "gridRefinement/octreeGridStructure.h"

namespace plb {

// Contains all the necessary information for a complete grid refined lattice.
template <typename T, int nDim>
class MultiLevelTensorField3D : public MultiLevel3D {
public:
    // The dynamics dyn must be right for the coarsest instantiated level.
    MultiLevelTensorField3D(const OctreeGridStructure &ogs_);
    MultiLevelTensorField3D(const OctreeGridStructure &ogs_, Box3D domain, plint domainLevel);

    MultiLevelTensorField3D(const OctreeGridStructure &ogs_, const Array<T, nDim> &iniVal);
    MultiLevelTensorField3D(
        const OctreeGridStructure &ogs_, Box3D domain, plint domainLevel,
        const Array<T, nDim> &iniVal);

    MultiLevelTensorField3D(MultiLevelTensorField3D<T, nDim> const &rhs);

    MultiLevelTensorField3D<T, nDim> &operator=(MultiLevelTensorField3D<T, nDim> const &rhs);
    void swap(MultiLevelTensorField3D<T, nDim> &rhs);
    ~MultiLevelTensorField3D();

    virtual MultiTensorField3D<T, nDim> const &getLevel(plint iL) const
    {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return *gridLevels[iL];
    }

    virtual MultiTensorField3D<T, nDim> &getLevel(plint iL)
    {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return *gridLevels[iL];
    }

    plint getNumLevels() const
    {
        return ogs.getNumLevels();
    }

    OctreeGridStructure const &getOgs() const
    {
        return ogs;
    }

private:
    OctreeGridStructure ogs;

private:
    std::vector<MultiTensorField3D<T, nDim> *> gridLevels;
};

// Contains all the necessary information for a complete grid refined lattice.
template <typename T, int nDim>
class MultiLevelTensorFieldForOutput3D : public MultiLevel3D {
public:
    // The dynamics dyn must be right for the coarsest instantiated level.
    MultiLevelTensorFieldForOutput3D(const OctreeGridStructure &ogs_, bool crop);
    MultiLevelTensorFieldForOutput3D(
        const OctreeGridStructure &ogs_, Box3D domain, plint domainLevel, bool crop);
    MultiLevelTensorFieldForOutput3D(MultiLevelTensorFieldForOutput3D<T, nDim> const &rhs);

    MultiLevelTensorFieldForOutput3D<T, nDim> &operator=(
        MultiLevelTensorFieldForOutput3D<T, nDim> const &rhs);
    void swap(MultiLevelTensorFieldForOutput3D<T, nDim> &rhs);
    ~MultiLevelTensorFieldForOutput3D();

    virtual MultiTensorField3D<T, nDim> const &getLevel(plint iL) const
    {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return *gridLevels[iL];
    }

    virtual MultiTensorField3D<T, nDim> &getLevel(plint iL)
    {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return *gridLevels[iL];
    }

    plint getNumLevels() const
    {
        return ogs.getNumLevels();
    }

    OctreeGridStructure const &getOgs() const
    {
        return ogs;
    }

private:
    OctreeGridStructure ogs;

private:
    std::vector<MultiTensorField3D<T, nDim> *> gridLevels;
};

}  // namespace plb

#endif  // MULTI_LEVEL_TENSOR_FIELD_3D_H
