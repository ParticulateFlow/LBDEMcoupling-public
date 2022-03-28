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

#ifndef INTERPOLATIONS_3D_HH
#define INTERPOLATIONS_3D_HH

#include <vector>

#include "core/globalDefs.h"
#include "core/util.h"
#include "finiteDifference/interpolations3D.h"

namespace plb {

/* ******** Function linearInterpolationCoefficients ********************* */

template <typename T>
void linearInterpolationCoefficients(
    AtomicBlock3D const &block, Array<T, 3> const &position, std::vector<Dot3D> &cellPos,
    std::vector<T> &weights)
{
    cellPos.resize(8);
    cellPos[0] = Dot3D((plint)position[0], (plint)position[1], (plint)position[2]);
    cellPos[1] = Dot3D((plint)position[0], (plint)position[1], (plint)(position[2] + (T)1.0));
    cellPos[2] = Dot3D((plint)position[0], (plint)(position[1] + (T)1.0), (plint)position[2]);
    cellPos[3] =
        Dot3D((plint)position[0], (plint)(position[1] + (T)1.0), (plint)(position[2] + (T)1.0));
    cellPos[4] = Dot3D((plint)(position[0] + (T)1.0), (plint)position[1], (plint)position[2]);
    cellPos[5] =
        Dot3D((plint)(position[0] + (T)1.0), (plint)position[1], (plint)(position[2] + (T)1.0));
    cellPos[6] =
        Dot3D((plint)(position[0] + (T)1.0), (plint)(position[1] + (T)1.0), (plint)position[2]);
    cellPos[7] = Dot3D(
        (plint)(position[0] + (T)1.0), (plint)(position[1] + (T)1.0),
        (plint)(position[2] + (T)1.0));

    T u = position[0] - (T)cellPos[0].x;
    T v = position[1] - (T)cellPos[0].y;
    T w = position[2] - (T)cellPos[0].z;

    weights.resize(8);
    weights[0] = (1. - u) * (1. - v) * (1. - w);
    weights[1] = (1. - u) * (1. - v) * (w);
    weights[2] = (1. - u) * (v) * (1. - w);
    weights[3] = (1. - u) * (v) * (w);
    weights[4] = (u) * (1. - v) * (1. - w);
    weights[5] = (u) * (1. - v) * (w);
    weights[6] = (u) * (v) * (1. - w);
    weights[7] = (u) * (v) * (w);

    // Convert cell position to local coordinates.
    for (plint iPos = 0; iPos < 8; ++iPos) {
        cellPos[iPos] -= block.getLocation();
    }
}

template <typename T, plint nDim>
Array<T, nDim> linearInterpolateTensorField(
    TensorField3D<T, nDim> &tensorField, Array<T, 3> const &position)
{
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(tensorField, position, pos, weights);
    Array<T, nDim> vector;
    vector.resetToZero();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        vector += weights[iCell] * tensorField.get(pos[iCell].x, pos[iCell].y, pos[iCell].z);
    }

    return vector;
}

template <typename T>
void linearInterpolateNtensorField(
    NTensorField3D<T> &tensorField, Array<T, 3> const &position, std::vector<T> &result)
{
    plint n = tensorField.getNdim();
    result.resize(n);
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(tensorField, position, pos, weights);
    for (plint i = 0; i < n; ++i) {
        result[i] = T();
    }
    for (plint i = 0; i < n; ++i) {
        for (plint iCell = 0; iCell < 8; ++iCell) {
            result[i] +=
                weights[iCell] * tensorField.get(pos[iCell].x, pos[iCell].y, pos[iCell].z)[i];
        }
    }
}

template <typename T, plint nDim>
Array<T, nDim> predictorCorrectorTensorField(
    TensorField3D<T, nDim> &tensorField, Array<T, 3> const &position, T scaling)
{
    Array<T, 3> position1(position);
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(tensorField, position1, pos, weights);
    Array<T, nDim> vector1;
    vector1.resetToZero();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        vector1 +=
            weights[iCell] * tensorField.get(pos[iCell].x, pos[iCell].y, pos[iCell].z) * scaling;
    }

    Array<T, 3> position2(position1 + vector1);
    linearInterpolationCoefficients(tensorField, position2, pos, weights);
    Array<T, nDim> vector2;
    vector2.resetToZero();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        vector2 +=
            weights[iCell] * tensorField.get(pos[iCell].x, pos[iCell].y, pos[iCell].z) * scaling;
    }

    return (vector1 + vector2) / (T)2;
}

template <typename T>
Array<T, 3> predictorCorrectorNTensorField(
    NTensorField3D<T> &tensorField, Array<T, 3> const &position, T scaling)
{
    PLB_PRECONDITION(tensorField.getNdim());
    Array<T, 3> position1(position);
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(tensorField, position1, pos, weights);
    Array<T, 3> vector1;
    vector1.resetToZero();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        T *data = tensorField.get(pos[iCell].x, pos[iCell].y, pos[iCell].z);
        vector1[0] += data[0] * weights[iCell] * scaling;
        vector1[1] += data[1] * weights[iCell] * scaling;
        vector1[2] += data[2] * weights[iCell] * scaling;
    }

    Array<T, 3> position2(position1 + vector1);
    linearInterpolationCoefficients(tensorField, position2, pos, weights);
    Array<T, 3> vector2;
    vector2.resetToZero();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        T *data = tensorField.get(pos[iCell].x, pos[iCell].y, pos[iCell].z);
        vector2[0] += data[0] * weights[iCell] * scaling;
        vector2[1] += data[1] * weights[iCell] * scaling;
        vector2[2] += data[2] * weights[iCell] * scaling;
    }

    return (vector1 + vector2) / (T)2;
}

template <typename T>
void predictorCorrectorRhoBarJ(
    NTensorField3D<T> &rhoBarJ, Array<T, 3> const &position, bool velIsJ, Array<T, 3> &j, T &rhoBar)
{
    PLB_ASSERT(rhoBarJ.getNdim() == 4);
    Array<T, 3> position1(position);
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(rhoBarJ, position1, pos, weights);
    Array<T, 3> j1;
    j1.resetToZero();
    T rhoBar1 = T();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        T const *data = rhoBarJ.get(pos[iCell].x, pos[iCell].y, pos[iCell].z);
        j1.add_from_cArray(data + 1, weights[iCell]);
        rhoBar1 += weights[iCell] * (*data);
    }

    Array<T, 3> position2(position1 + j1);
    linearInterpolationCoefficients(rhoBarJ, position2, pos, weights);
    Array<T, 3> j2;
    j2.resetToZero();
    T rhoBar2 = T();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        T const *data = rhoBarJ.get(pos[iCell].x, pos[iCell].y, pos[iCell].z);
        j2.add_from_cArray(data + 1, weights[iCell]);
        rhoBar2 += weights[iCell] * (*data);
    }

    j = (j1 + j2) / (T)2;
    rhoBar = (rhoBar1 + rhoBar2) / (T)2;
}

template <typename T>
T linearInterpolateScalarField(ScalarField3D<T> &scalarField, Array<T, 3> const &position)
{
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(scalarField, position, pos, weights);
    T scalar = T();
    for (plint iCell = 0; iCell < 8; ++iCell) {
        scalar += weights[iCell] * scalarField.get(pos[iCell].x, pos[iCell].y, pos[iCell].z);
    }

    return scalar;
}

}  // namespace plb

#endif  // INTERPOLATIONS_3D_HH
