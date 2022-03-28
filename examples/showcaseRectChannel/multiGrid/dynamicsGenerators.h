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
 * Generator of Dynamics vectors for multigrid
 */
#ifndef DYNAMICS_GENERATORS_H
#define DYNAMICS_GENERATORS_H

#include "basicDynamics/isoThermalDynamics.h"
#include "multiGrid/multiGridParameterManager.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
std::vector<Dynamics<T, Descriptor> *> generateBGKdynamics(
    RefinementParameters<T> *refinementParameters)
{
    plint numLevel = refinementParameters->getNumLevels();

    std::vector<Dynamics<T, Descriptor> *> result(numLevel);
    for (plint iLevel = 0; iLevel < numLevel; ++iLevel) {
        T omega = refinementParameters->getParameters(iLevel).getOmega();
        result[iLevel] = new BGKdynamics<T, Descriptor>(omega);
    }

    return result;
}

template <typename T, template <typename U> class Descriptor>
std::vector<Dynamics<T, Descriptor> *> generateBounceBackDynamics(
    RefinementParameters<T> *refinementParameters)
{
    plint numLevel = refinementParameters->getNumLevels();

    std::vector<Dynamics<T, Descriptor> *> result(numLevel);
    for (plint iLevel = 0; iLevel < numLevel; ++iLevel) {
        result[iLevel] = new BounceBack<T, Descriptor>(0.0);
    }

    return result;
}

// TODO add generators for all dynamics that will be used

}  // namespace plb

#endif  // DYNAMICS_GENERATORS_H
