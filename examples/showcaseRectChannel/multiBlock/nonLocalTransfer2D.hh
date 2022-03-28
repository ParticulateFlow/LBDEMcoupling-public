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
 * Geometry specifications for 2D multiblocks -- header file.
 */

#ifndef NON_LOCAL_TRANSFER_2D_HH
#define NON_LOCAL_TRANSFER_2D_HH

#include <vector>

#include "core/globalDefs.h"
#include "multiBlock/nonLocalTransfer2D.h"

namespace plb {

template <typename T>
void copyNonLocal(MultiScalarField2D<T> const &from, MultiScalarField2D<T> &to, Box2D const &domain)
{
    // copyNonLocal_generic(from, to, domain, modif::staticVariables);
    copy(from, domain, to, domain);
}

template <typename T>
void copy(
    MultiScalarField2D<T> const &from, Box2D const &fromDomain, MultiScalarField2D<T> &to,
    Box2D const &toDomain)
{
    copy_generic(from, fromDomain, to, toDomain, modif::staticVariables);
}

template <typename T>
void copyNonLocal(
    MultiNTensorField2D<T> const &from, MultiNTensorField2D<T> &to, Box2D const &domain)
{
    // copyNonLocal_generic(from, to, domain, modif::staticVariables);
    copy(from, domain, to, domain);
}

template <typename T>
void copy(
    MultiNTensorField2D<T> const &from, Box2D const &fromDomain, MultiNTensorField2D<T> &to,
    Box2D const &toDomain)
{
    copy_generic(from, fromDomain, to, toDomain, modif::staticVariables);
}

template <typename T, int nDim>
void copyNonLocal(
    MultiTensorField2D<T, nDim> const &from, MultiTensorField2D<T, nDim> &to, Box2D const &domain)
{
    // copyNonLocal_generic(from, to, domain, modif::staticVariables);
    copy(from, domain, to, domain);
}

template <typename T, int nDim>
void copy(
    MultiTensorField2D<T, nDim> const &from, Box2D const &fromDomain,
    MultiTensorField2D<T, nDim> &to, Box2D const &toDomain)
{
    copy_generic(from, fromDomain, to, toDomain, modif::staticVariables);
}

template <typename T, template <typename U> class Descriptor>
void copyNonLocal(
    MultiBlockLattice2D<T, Descriptor> const &from, MultiBlockLattice2D<T, Descriptor> &to,
    Box2D const &domain, modif::ModifT whichContent)
{
    // copyNonLocal_generic(from, to, domain, whichContent);
    copy(from, domain, to, domain, whichContent);
}

template <typename T, template <typename U> class Descriptor>
void copy(
    MultiBlockLattice2D<T, Descriptor> const &from, Box2D const &fromDomain,
    MultiBlockLattice2D<T, Descriptor> &to, Box2D const &toDomain, modif::ModifT whichContent)
{
    copy_generic(from, fromDomain, to, toDomain, whichContent);
}

template <typename T, template <typename U> class Descriptor>
void copyPopulations(
    MultiBlockLattice2D<T, Descriptor> const &from, Box2D const &fromDomain,
    MultiBlockLattice2D<T, Descriptor> &to, Box2D const &toDomain)
{
    copy(from, fromDomain, to, toDomain, modif::staticVariables);
}

template <typename T, template <typename U> class Descriptor>
void copyAll(
    MultiBlockLattice2D<T, Descriptor> const &from, Box2D const &fromDomain,
    MultiBlockLattice2D<T, Descriptor> &to, Box2D const &toDomain)
{
    copy(from, fromDomain, to, toDomain, modif::allVariables);
}

template <typename T, template <typename U> class Descriptor>
void copyRegenerate(
    MultiBlockLattice2D<T, Descriptor> const &from, Box2D const &fromDomain,
    MultiBlockLattice2D<T, Descriptor> &to, Box2D const &toDomain)
{
    copy(from, fromDomain, to, toDomain, modif::dataStructure);
}

}  // namespace plb

#endif  // NON_LOCAL_TRANSFER_2D_HH