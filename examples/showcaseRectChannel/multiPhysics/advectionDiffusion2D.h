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

#ifndef ADVECTION_DIFFUSION_2D_H
#define ADVECTION_DIFFUSION_2D_H

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/block2D.h"
#include "core/globalDefs.h"

namespace plb {

/**
 * Multiphysics class for one-way coupling between Navier-Stokes and
 * advection-diffusion equations: the fluid velocity is copied
 * to the advection-diffusion field, which is advected passively.
 */
template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
class LatticeToPassiveAdvDiff2D :
    public BoxProcessingFunctional2D_LL<T, FluidDescriptor, T, ScalarDescriptor> {
public:
    LatticeToPassiveAdvDiff2D(T scaling_ = 1.);
    virtual void process(
        Box2D domain, BlockLattice2D<T, FluidDescriptor> &fluid,
        BlockLattice2D<T, ScalarDescriptor> &scalar);
    virtual LatticeToPassiveAdvDiff2D<T, FluidDescriptor, ScalarDescriptor> *clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T scaling;
};

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void latticeToPassiveAdvDiff(
    MultiBlockLattice2D<T, FluidDescriptor> &fluid,
    MultiBlockLattice2D<T, ScalarDescriptor> &scalar, Box2D domain);

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_2D_H
