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
 * 3D couplings -- generic implementation
 */
#ifndef COUPLING_3D_HH
#define COUPLING_3D_HH

#include "core/globalDefs.h"
#include "multiBlock/coupling3D.h"

namespace plb {

/* ********************  FullDomainCollideAndStreamAction3D ************************** */

template <typename T, template <typename U> class Descriptor>
FullDomainCollideAndStreamAction3D<T, Descriptor>::FullDomainCollideAndStreamAction3D(
    plint blockId_) :
    blockId(blockId_)
{ }

template <typename T, template <typename U> class Descriptor>
FullDomainCollideAndStreamAction3D<T, Descriptor>
    *FullDomainCollideAndStreamAction3D<T, Descriptor>::clone() const
{
    return new FullDomainCollideAndStreamAction3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void FullDomainCollideAndStreamAction3D<T, Descriptor>::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock3D *block = multiBlockRegistration3D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    MultiBlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice3D<T, Descriptor> *>(block);
    PLB_ASSERT(lattice);
    lattice->externalCollideAndStream();
}

/* ********************  CollideAndStreamAction3D ************************** */

template <typename T, template <typename U> class Descriptor>
CollideAndStreamAction3D<T, Descriptor>::CollideAndStreamAction3D(plint blockId_, Box3D domain_) :
    blockId(blockId_), domain(domain_)
{ }

template <typename T, template <typename U> class Descriptor>
CollideAndStreamAction3D<T, Descriptor> *CollideAndStreamAction3D<T, Descriptor>::clone() const
{
    return new CollideAndStreamAction3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CollideAndStreamAction3D<T, Descriptor>::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock3D *block = multiBlockRegistration3D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    MultiBlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice3D<T, Descriptor> *>(block);
    PLB_ASSERT(lattice);
    lattice->collideAndStream(domain);
}

/* ********************  FullDomainStreamAction3D ************************** */

template <typename T, template <typename U> class Descriptor>
FullDomainStreamAction3D<T, Descriptor>::FullDomainStreamAction3D(plint blockId_) :
    blockId(blockId_)
{ }

template <typename T, template <typename U> class Descriptor>
FullDomainStreamAction3D<T, Descriptor> *FullDomainStreamAction3D<T, Descriptor>::clone() const
{
    return new FullDomainStreamAction3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void FullDomainStreamAction3D<T, Descriptor>::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock3D *block = multiBlockRegistration3D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    MultiBlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice3D<T, Descriptor> *>(block);
    PLB_ASSERT(lattice);
    lattice->externalStream();
}

/* ********************  StreamAction3D ************************** */

template <typename T, template <typename U> class Descriptor>
StreamAction3D<T, Descriptor>::StreamAction3D(plint blockId_, Box3D domain_) :
    blockId(blockId_), domain(domain_)
{ }

template <typename T, template <typename U> class Descriptor>
StreamAction3D<T, Descriptor> *StreamAction3D<T, Descriptor>::clone() const
{
    return new StreamAction3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void StreamAction3D<T, Descriptor>::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock3D *block = multiBlockRegistration3D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    MultiBlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice3D<T, Descriptor> *>(block);
    PLB_ASSERT(lattice);
    lattice->stream(domain);
}

/* ********************  IncrementTimeAction3D ************************** */

template <typename T, template <typename U> class Descriptor>
IncrementTimeAction3D<T, Descriptor>::IncrementTimeAction3D(plint blockId_) : blockId(blockId_)
{ }

template <typename T, template <typename U> class Descriptor>
IncrementTimeAction3D<T, Descriptor> *IncrementTimeAction3D<T, Descriptor>::clone() const
{
    return new IncrementTimeAction3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void IncrementTimeAction3D<T, Descriptor>::execute(std::vector<id_t> &allMultiBlocks)
{
    PLB_ASSERT(blockId < (plint)allMultiBlocks.size());
    MultiBlock3D *block = multiBlockRegistration3D().find(allMultiBlocks[blockId]);
    PLB_ASSERT(block);
    MultiBlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<MultiBlockLattice3D<T, Descriptor> *>(block);
    PLB_ASSERT(lattice);
    lattice->incrementTime();
}

/* ********************  ExecuteFunctionAction ************************** */

template <class Function>
ExecuteFunctionAction<Function>::ExecuteFunctionAction(Function f_) : f(f_)
{ }

template <class Function>
ExecuteFunctionAction<Function> *ExecuteFunctionAction<Function>::clone() const
{
    return new ExecuteFunctionAction<Function>(*this);
}

template <class Function>
void ExecuteFunctionAction<Function>::execute(std::vector<id_t> &allMultiBlocks)
{
    f();
}

/* ********************  Actions3D ************************** */

template <typename T, template <typename U> class Descriptor>
plint Actions3D::addCollideAndStream(plint blockNum, Box3D domain)
{
    return addAction(new CollideAndStreamAction3D<T, Descriptor>(blockNum, domain));
}

template <typename T, template <typename U> class Descriptor>
plint Actions3D::addCollideAndStream(plint blockNum)
{
    return addAction(new FullDomainCollideAndStreamAction3D<T, Descriptor>(blockNum));
}

template <typename T, template <typename U> class Descriptor>
plint Actions3D::addStream(plint blockNum, Box3D domain)
{
    return addAction(new StreamAction3D<T, Descriptor>(blockNum, domain));
}

template <typename T, template <typename U> class Descriptor>
plint Actions3D::addStream(plint blockNum)
{
    return addAction(new FullDomainStreamAction3D<T, Descriptor>(blockNum));
}

template <typename T, template <typename U> class Descriptor>
plint Actions3D::addIncrementTime(plint blockNum)
{
    return addAction(new IncrementTimeAction3D<T, Descriptor>(blockNum));
}

template <class Function>
plint Actions3D::addExecuteFunction(Function f)
{
    return addAction(new ExecuteFunctionAction<Function>(f));
}

}  // namespace plb

#endif  // COUPLING_3D_HH
