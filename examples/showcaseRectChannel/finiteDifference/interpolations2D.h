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

#ifndef INTERPOLATIONS_2D_H
#define INTERPOLATIONS_2D_H

#include <vector>

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/dataField2D.h"
#include "core/geometry2D.h"
#include "core/globalDefs.h"

namespace plb {

/// Helper function: linear interpolation within one cell.
template <typename T>
void linearInterpolationCoefficients(AtomicBlock2D const &block, Array<T, 2> const &position);

template <typename T, plint nDim>
Array<T, nDim> linearInterpolateTensorField(
    TensorField2D<T, nDim> &tensorField, Array<T, 2> const &position);

template <typename T, plint nDim>
Array<T, nDim> predictorCorrectorTensorField(
    TensorField2D<T, nDim> &tensorField, Array<T, 2> const &position, T scaling);

template <typename T>
Array<T, 2> predictorCorrectorNTensorField(
    NTensorField2D<T> &tensorField, Array<T, 2> const &position, T scaling);

template <typename T>
void predictorCorrectorRhoBarJ(
    NTensorField2D<T> &rhoBarJ, Array<T, 2> const &position, bool velIsJ, Array<T, 2> &j,
    T &rhoBar);

template <typename T>
T linearInterpolateScalarField(ScalarField2D<T> &scalarField, Array<T, 2> const &position);

}  // namespace plb

#endif  // INTERPOLATIONS_2D_H
