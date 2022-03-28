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
 * Helper functions for domain initialization -- header file.
 */
#ifndef THERMAL_DATA_ANALYSIS_3D_HH
#define THERMAL_DATA_ANALYSIS_3D_HH

#include <cmath>

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "core/blockStatistics.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "multiPhysics/thermalDataAnalysis3D.h"

namespace plb {

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
T computeNusseltNumber(
    BlockLattice3D<T, FluidDescriptor> &fluid,
    BlockLattice3D<T, TemperatureDescriptor> &temperature, Box3D domain, int direction, T deltaX,
    T kappa, T deltaTemperature)
{
    BoxSumVelocityTemperatureFunctional3D<T, FluidDescriptor, TemperatureDescriptor> functional(
        direction);
    applyProcessingFunctional(functional, domain, fluid, temperature);
    return (T)1
           + functional.getSumVelocityTemperature()
                 / ((T)domain.nCells() * deltaX * deltaTemperature * kappa);
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
T computeNusseltNumber(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, TemperatureDescriptor> &temperature, Box3D domain, int direction,
    T deltaX, T kappa, T deltaTemperature)
{
    BoxSumVelocityTemperatureFunctional3D<T, FluidDescriptor, TemperatureDescriptor> functional(
        direction);
    applyProcessingFunctional(functional, domain, fluid, temperature);
    return (T)1
           + functional.getSumVelocityTemperature()
                 / ((T)domain.nCells() * deltaX * deltaTemperature * kappa);
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
BoxSumVelocityTemperatureFunctional3D<T, FluidDescriptor, TemperatureDescriptor>::
    BoxSumVelocityTemperatureFunctional3D(int direction_) :
    direction(direction_), sumVelocityTemperatureId(this->getStatistics().subscribeSum())
{ }

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
void BoxSumVelocityTemperatureFunctional3D<T, FluidDescriptor, TemperatureDescriptor>::process(
    Box3D domain, BlockLattice3D<T, FluidDescriptor> &fluid,
    BlockLattice3D<T, TemperatureDescriptor> &temperature)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, FluidDescriptor<T>::d> velocity;
                fluid.get(iX, iY, iZ).computeVelocity(velocity);
                T localTemperature = temperature.get(iX, iY, iZ).computeDensity();

                T velocityTemperature = localTemperature * velocity[direction];
                statistics.gatherSum(sumVelocityTemperatureId, velocityTemperature);
            }
        }
    }
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
BoxSumVelocityTemperatureFunctional3D<T, FluidDescriptor, TemperatureDescriptor>
    *BoxSumVelocityTemperatureFunctional3D<T, FluidDescriptor, TemperatureDescriptor>::clone() const
{
    return new BoxSumVelocityTemperatureFunctional3D<T, FluidDescriptor, TemperatureDescriptor>(
        *this);
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
T BoxSumVelocityTemperatureFunctional3D<
    T, FluidDescriptor, TemperatureDescriptor>::getSumVelocityTemperature() const
{
    return this->getStatistics().getSum(sumVelocityTemperatureId);
}

}  // namespace plb

#endif  // THERMAL_DATA_ANALYSIS_3D_HH
