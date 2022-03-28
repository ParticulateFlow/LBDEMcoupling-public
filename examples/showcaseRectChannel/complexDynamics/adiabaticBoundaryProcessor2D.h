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
 * Neumann boundary conditions for temperature -- header file.
 */
#ifndef ADIABATIC_BOUNDARY_PROCESSOR_2D_H
#define ADIABATIC_BOUNDARY_PROCESSOR_2D_H

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "core/blockLatticeBase2D.h"
#include "core/globalDefs.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class FlatAdiabaticBoundaryFunctional2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    virtual plint extent() const
    {
        return 2;
    }
    virtual plint extent(int whichDirection) const
    {
        return 2;
    }
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual FlatAdiabaticBoundaryFunctional2D<T, Descriptor, direction, orientation> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // ADIABATIC_BOUNDARY_PROCESSOR_2D_H
