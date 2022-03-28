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

#ifndef BOUSSINESQ_THERMAL_PROCESSOR_2D_H
#define BOUSSINESQ_THERMAL_PROCESSOR_2D_H

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/block2D.h"
#include "core/globalDefs.h"

namespace plb {

/**
 * Class for the coupling between a Navier-Stokes (NS) lattice and an
 * Advection-Diffusion (AD) lattice in the boussinesq approximation.
 */
template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
class BoussinesqThermalProcessor2D :
    public BoxProcessingFunctional2D_LL<T, FluidDescriptor, T, TemperatureDescriptor> {
public:
    BoussinesqThermalProcessor2D(
        T gravity_, T T0_, T deltaTemp_, Array<T, FluidDescriptor<T>::d> dir_);

    virtual void process(
        Box2D domain, BlockLattice2D<T, FluidDescriptor> &fluid,
        BlockLattice2D<T, TemperatureDescriptor> &temperature);
    virtual BoussinesqThermalProcessor2D<T, FluidDescriptor, TemperatureDescriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
        modified[1] = modif::staticVariables;
    }

private:
    T gravity, T0, deltaTemp;
    Array<T, FluidDescriptor<T>::d> dir;
};

/**
 * Class for the coupling between a Navier-Stokes (NS) lattice and an
 * Advection-Diffusion (AD) lattice in the boussinesq approximation.
 */
template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class TemperatureDescriptor>
class CompleteBoussinesqThermalProcessor2D :
    public BoxProcessingFunctional2D_LL<T, FluidDescriptor, T, TemperatureDescriptor> {
public:
    CompleteBoussinesqThermalProcessor2D(
        T gravity_, T T0_, T deltaTemp_, Array<T, FluidDescriptor<T>::d> dir_);

    virtual void process(
        Box2D domain, BlockLattice2D<T, FluidDescriptor> &fluid,
        BlockLattice2D<T, TemperatureDescriptor> &temperature);
    virtual CompleteBoussinesqThermalProcessor2D<T, FluidDescriptor, TemperatureDescriptor> *clone()
        const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
        modified[1] = modif::staticVariables;
    }

private:
    T gravity, T0, deltaTemp;
    Array<T, FluidDescriptor<T>::d> dir;
};

}  // namespace plb

#endif  // BOUSSINESQ_THERMAL_PROCESSOR_2D_H
