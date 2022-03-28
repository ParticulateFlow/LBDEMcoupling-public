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

#ifndef INTERPOLATION_HELPER_HH
#define INTERPOLATION_HELPER_HH

#include "interpolationHelper.h"

namespace plb {

template <typename T>
std::vector<T> asymetricCubicInterpolation(T f[4][3])
{
    std::vector<T> result(4);

    result[0] = 3. / 8. * f[1][0] + 3. / 4. * f[1][1] - 1. / 8. * f[1][2];
    result[1] = 9. / 16. * (f[1][0] + f[2][0]) - 1. / 16. * (f[0][0] + f[3][0]);
    result[2] = -(3. / 128.) * f[0][0] - (3. / 64.) * f[0][1] + (1. / 128.) * f[0][2]
                + (27. / 128.) * f[1][0] + (27. / 64.) * f[1][1] - (9. / 128.) * f[1][2]
                + (27. / 128.) * f[2][0] + (27. / 64.) * f[2][1] - (9. / 128.) * f[2][2]
                - (3. / 128.) * f[3][0] - (3. / 64.) * f[3][1] + (1. / 128.) * f[3][2];
    result[3] = 3. / 8. * f[2][0] + 3. / 4. * f[2][1] - 1. / 8. * f[2][2];

    return result;
}

template <typename T>
std::vector<T> symetricCubicInterpolation(T f[4][4])
{
    std::vector<T> result(5);

    result[0] = 9. / 16. * (f[1][1] + f[1][2]) - 1. / 16. * (f[1][0] + f[1][3]);
    result[1] = 9. / 16. * (f[1][1] + f[2][1]) - 1. / 16. * (f[0][1] + f[3][1]);

    result[2] = (1. / 256.) * f[0][0] - (9. / 256.) * f[0][1] - (9. / 256.) * f[0][2]
                + (1. / 256.) * f[0][3] - (9. / 256.) * f[1][0] + (81. / 256.) * f[1][1]
                + (81. / 256.) * f[1][2] - (9. / 256.) * f[1][3] - (9. / 256.) * f[2][0]
                + (81. / 256.) * f[2][1] + (81. / 256.) * f[2][2] - (9. / 256.) * f[2][3]
                + (1. / 256.) * f[3][0] - (9. / 256.) * f[3][1] - (9. / 256.) * f[3][2]
                + (1. / 256.) * f[3][3];

    result[3] = 9. / 16. * (f[1][2] + f[2][2]) - 1. / 16. * (f[0][2] + f[3][2]);
    result[4] = 9. / 16. * (f[2][1] + f[2][2]) - 1. / 16. * (f[2][0] + f[2][3]);

    return result;
}

template <typename T>
std::vector<T> cornerInterpolation(T f[3][3])
{
    std::vector<T> result(3);
    result[0] = 3. / 8. * f[0][0] + 3. / 4. * f[0][1] - 1. / 8. * f[0][2];
    result[1] = 3. / 8. * f[0][0] + 3. / 4. * f[1][0] - 1. / 8. * f[2][0];
    result[2] = (9. / 64.) * f[0][0] + (9. / 32.) * f[0][1] - (3. / 64.) * f[0][2]
                + (9. / 32.) * f[1][0] + (9. / 16.) * f[1][1] - (3. / 32.) * f[1][2]
                - (3. / 64.) * f[2][0] - (3. / 32.) * f[2][1] + (1. / 64.) * f[2][2];

    return result;
}

template <typename T>
std::vector<T> helperCornerInterpolation(T f[3][3])
{
    std::vector<T> result(5);
    result[0] = 3. / 8. * f[0][0] + 3. / 4. * f[0][1] - 1. / 8. * f[0][2];
    result[1] = 3. / 8. * f[0][0] + 3. / 4. * f[1][0] - 1. / 8. * f[2][0];
    result[2] = (9. / 64.) * f[0][0] + (9. / 32.) * f[0][1] - (3. / 64.) * f[0][2]
                + (9. / 32.) * f[1][0] + (9. / 16.) * f[1][1] - (3. / 32.) * f[1][2]
                - (3. / 64.) * f[2][0] - (3. / 32.) * f[2][1] + (1. / 64.) * f[2][2];
    result[3] = 3. / 8. * f[0][1] + 3. / 4. * f[1][1] - 1. / 8. * f[2][1];
    result[4] = 3. / 8. * f[1][0] + 3. / 4. * f[1][1] - 1. / 8. * f[1][2];

    return result;
}

/// Function to copy the populations to a cell of the fine grid
template <typename T, template <typename U> class Descriptor>
void copyPopulations(std::vector<T> &decomposedValues, Cell<T, Descriptor> &cell)
{
    plint whichTime = 1;
    dynamic_cast<FineGridBoundaryDynamics<T, Descriptor> &>(cell.getDynamics())
        .getDecomposedValues(whichTime)
        .assign(decomposedValues.begin(), decomposedValues.end());
}

}  // namespace plb

#endif  // INTERPOLATION_HELPER_HH
