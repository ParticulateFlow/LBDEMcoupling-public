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
#ifndef THERMAL_DATA_ANALYSIS_3D_H
#define THERMAL_DATA_ANALYSIS_3D_H

#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/reductiveDataProcessorWrapper3D.h"
#include "core/blockLatticeBase3D.h"
#include "core/dataFieldBase3D.h"
#include "core/globalDefs.h"

namespace plb {

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
T computeNusseltNumber(
    BlockLattice3D<T, FluidDescriptor> &fluid,
    BlockLattice3D<T, TemperatureDescriptor> &temperature, Box3D domain, int direction, T deltaX,
    T kappa, T deltaTemperature);

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
T computeNusseltNumber(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, TemperatureDescriptor> &temperature, Box3D domain, int direction,
    T deltaX, T kappa, T deltaTemperature);

/** This class computes the product between the component "direction" of the
 *  velocity (computed from the fluid-lattice) and the temperature
 *  computed from the advection-diffusion-lattice), and returns its sum
 *  over the whole domain.
 */
template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
class BoxSumVelocityTemperatureFunctional3D :
    public ReductiveBoxProcessingFunctional3D_LL<T, FluidDescriptor, T, TemperatureDescriptor> {
public:
    BoxSumVelocityTemperatureFunctional3D(int direction_);
    virtual void process(
        Box3D domain, BlockLattice3D<T, FluidDescriptor> &fluid,
        BlockLattice3D<T, TemperatureDescriptor> &temperature);
    virtual BoxSumVelocityTemperatureFunctional3D<T, FluidDescriptor, TemperatureDescriptor>
        *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getSumVelocityTemperature() const;

private:
    int direction;
    plint sumVelocityTemperatureId;
};

}  // namespace plb

#endif  // THERMAL_DATA_ANALYSIS_3D_H
