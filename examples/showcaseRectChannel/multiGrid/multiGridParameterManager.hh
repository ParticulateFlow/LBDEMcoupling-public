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

#ifndef MULTI_GRID_PARAMETER_MANAGER_HH
#define MULTI_GRID_PARAMETER_MANAGER_HH

#include "multiGrid/multiGridParameterManager.h"

namespace plb {

/// Interface for a wrapper of refinement parameters (each level posseses specific parameters)
template <typename T>
RefinementParameters<T>::RefinementParameters(
    plint levelNumber_, plint referenceLevel_, IncomprFlowParam<T> parameters_) :
    levelNumber(levelNumber_), referenceLevel(referenceLevel_), originalParameters(parameters_)
{ }

template <typename T>
RefinementParameters<T>::RefinementParameters(RefinementParameters<T> const &rhs) :
    levelNumber(rhs.levelNumber),
    referenceLevel(rhs.referenceLevel),
    originalParameters(rhs.parameters_)
{ }

template <typename T>
IncomprFlowParam<T> const &RefinementParameters<T>::getParameters(plint lvl) const
{
    PLB_PRECONDITION(lvl >= 0 && lvl < (plint)parameters.size());
    return parameters[lvl];
}
template <typename T>
plint RefinementParameters<T>::getReferenceLevel()
{
    return referenceLevel;
}

template <typename T>
plint RefinementParameters<T>::getNumLevels()
{
    return levelNumber;
}

template <typename T>
IncomprFlowParam<T> const &RefinementParameters<T>::operator[](plint lvl) const
{
    return getParameters(lvl);
}

template <typename T>
void RefinementParameters<T>::putParameter(IncomprFlowParam<T> &newParam)
{
    parameters.push_back(newParam);
}

template <typename T>
IncomprFlowParam<T> &RefinementParameters<T>::getOriginalParameters()
{
    return originalParameters;
}

template <typename T>
ConvectiveRefinementParameters<T>::ConvectiveRefinementParameters(
    plint levels, plint reference, IncomprFlowParam<T> params) :
    RefinementParameters<T>(levels, reference, params)
{
    createParameters();
}

template <typename T>
void ConvectiveRefinementParameters<T>::createParameters()
{
    IncomprFlowParam<T> referenceParameters = RefinementParameters<T>::getOriginalParameters();
    if (RefinementParameters<T>::getReferenceLevel() == 0) {
        plint resolutionAtThisLevel = referenceParameters.getResolution();
        for (int iLvl = 0; iLvl < RefinementParameters<T>::getNumLevels(); ++iLvl) {
            // everything is the same but the resolution
            IncomprFlowParam<T> thisLevelParameters(
                referenceParameters.getLatticeU(), referenceParameters.getRe(),
                resolutionAtThisLevel, referenceParameters.getLx(), referenceParameters.getLy(),
                referenceParameters.getLz());
            this->putParameter(thisLevelParameters);
            resolutionAtThisLevel = 2 * resolutionAtThisLevel;
        }
    } else {
        // TODO implement creation of parameters if the referenceLevel != 0
        PLB_ASSERT(false);
    }
}

template <typename T>
DiffusiveRefinementParameters<T>::DiffusiveRefinementParameters(
    plint levels, plint reference, IncomprFlowParam<T> params) :
    RefinementParameters<T>(levels, reference)
{
    createParameters();
}

template <typename T>
void DiffusiveRefinementParameters<T>::createParameters()
{
    IncomprFlowParam<T> referenceParameters = RefinementParameters<T>::getOriginalParameters();
    if (RefinementParameters<T>::getReferenceLevel() == 0) {
        plint resolutionAtThisLevel = referenceParameters.getResolution();
        plint uMaxRef = referenceParameters.getLatticeU();
        plint Nref = referenceParameters.getResolution();
        for (int iLvl = 0; iLvl < RefinementParameters<T>::getNumLevels(); ++iLvl) {
            T uMax = uMaxRef * (T)resolutionAtThisLevel / (T)Nref;

            IncomprFlowParam<T> thisLevelParameters(
                uMax, referenceParameters.getRe(), resolutionAtThisLevel,
                referenceParameters.getLx(), referenceParameters.getLy(),
                referenceParameters.getLz());

            this->putParameter(thisLevelParameters);
            resolutionAtThisLevel = 2 * resolutionAtThisLevel;
        }
    } else {
        // TODO implement creation of parameters if the referenceLevel != 0
        PLB_ASSERT(false);
    }
}

}  // namespace plb

#endif  // MULTI_GRID_PARAMETER_MANAGER_H
