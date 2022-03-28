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
#ifndef META_STUFF_FUNCTIONAL_2D_H
#define META_STUFF_FUNCTIONAL_2D_H

#include <map>

#include "atomicBlock/atomicContainerBlock2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "core/globalDefs.h"
#include "dataProcessors/metaStuffHelper.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class StoreDynamicsFunctional2D : public PlainReductiveBoxProcessingFunctional2D {
public:
    StoreDynamicsFunctional2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual StoreDynamicsFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    pluint getMaxChainLength() const;

private:
    plint maxChainLengthId;
};

class IterateDynamicsFunctional2D : public PlainReductiveBoxProcessingFunctional2D {
public:
    IterateDynamicsFunctional2D(std::vector<int> previousMaximum_);
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual IterateDynamicsFunctional2D *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    std::vector<int> getNextMaximum() const;

private:
    std::vector<int> previousMaximum;
    std::vector<plint> maxIds;
};

template <typename T, template <typename U> class Descriptor>
class ExtractDynamicsChainFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, int> {
public:
    typedef std::map<std::vector<int>, int, VectorIsLess> DMap;

public:
    ExtractDynamicsChainFunctional2D(DMap const &dynamicsMap_, pluint maxChainSize_);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<int> &mask);
    virtual ExtractDynamicsChainFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    DMap dynamicsMap;
    pluint maxChainSize;
};

template <typename T, template <typename U> class Descriptor>
class ExtractTopMostDynamicsFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, int> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<int> &mask);
    virtual ExtractTopMostDynamicsFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class ExtractBottomMostDynamicsFunctional2D :
    public BoxProcessingFunctional2D_LS<T, Descriptor, int> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<int> &mask);
    virtual ExtractBottomMostDynamicsFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class AssignEntireCellFunctional2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor, T, Descriptor> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &sourceLattice,
        BlockLattice2D<T, Descriptor> &destinationLattice);
    virtual AssignEntireCellFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

class AllFlagsTrueFunctional2D : public PlainReductiveBoxProcessingFunctional2D {
public:
    AllFlagsTrueFunctional2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual AllFlagsTrueFunctional2D *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    bool allTrue() const;

private:
    plint numFalseId;
};

}  // namespace plb

#endif  // META_STUFF_FUNCTIONAL_2D_H
