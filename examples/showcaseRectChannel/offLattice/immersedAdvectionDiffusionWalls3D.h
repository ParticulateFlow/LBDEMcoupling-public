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

#ifndef IMMERSED_ADVECTION_DIFFUSION_WALLS_3D_H
#define IMMERSED_ADVECTION_DIFFUSION_WALLS_3D_H

#include <vector>

#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "immersedWalls3D.h"

namespace plb {

/* ******** InamuroAdvectionDiffusionIteration3D ************************************ */

template <typename T, class ScalarFunction>
class InamuroAdvectionDiffusionIteration3D : public BoxProcessingFunctional3D {
public:
    InamuroAdvectionDiffusionIteration3D(ScalarFunction scalarFunction_, T tau_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InamuroAdvectionDiffusionIteration3D<T, ScalarFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    ScalarFunction scalarFunction;
    T tau;
};

template <typename T, class ScalarFunction>
void inamuroAdvectionDiffusionIteration(
    ScalarFunction scalarFunction, MultiScalarField3D<T> &rhoBar, MultiContainerBlock3D &container,
    T tau)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&rhoBar);
    args.push_back(&container);
    applyProcessingFunctional(
        new InamuroAdvectionDiffusionIteration3D<T, ScalarFunction>(scalarFunction, tau),
        rhoBar.getBoundingBox(), args);
}

/* ******** IndexedInamuroAdvectionDiffusionIteration3D ************************************ */

// This is the same as InamuroAdvectionDiffusionIteration3D, with the difference that
// the ScalarFunction accepts as argument a global vertex index instead of
// a 3D position in space.
template <typename T, class ScalarFunction>
class IndexedInamuroAdvectionDiffusionIteration3D : public BoxProcessingFunctional3D {
public:
    IndexedInamuroAdvectionDiffusionIteration3D(ScalarFunction scalarFunction_, T tau_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual IndexedInamuroAdvectionDiffusionIteration3D<T, ScalarFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    ScalarFunction scalarFunction;
    T tau;
};

template <typename T, class ScalarFunction>
void indexedInamuroAdvectionDiffusionIteration(
    ScalarFunction scalarFunction, MultiScalarField3D<T> &rhoBar, MultiContainerBlock3D &container,
    T tau)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&rhoBar);
    args.push_back(&container);
    applyProcessingFunctional(
        new IndexedInamuroAdvectionDiffusionIteration3D<T, ScalarFunction>(scalarFunction, tau),
        rhoBar.getBoundingBox(), args);
}

/* ******** ConstScalarInamuroAdvectionDiffusionIteration3D ************************************ */

template <typename T>
class ConstScalarInamuroAdvectionDiffusionIteration3D : public BoxProcessingFunctional3D {
public:
    ConstScalarInamuroAdvectionDiffusionIteration3D(T wallScalar_, T tau_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ConstScalarInamuroAdvectionDiffusionIteration3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T wallScalar;
    T tau;
};

template <typename T>
void constScalarInamuroAdvectionDiffusionIteration(
    T wallScalar, MultiScalarField3D<T> &rhoBar, MultiContainerBlock3D &container, T tau)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&rhoBar);
    args.push_back(&container);
    applyProcessingFunctional(
        new ConstScalarInamuroAdvectionDiffusionIteration3D<T>(wallScalar, tau),
        rhoBar.getBoundingBox(), args);
}

/* ******** InstantiateImmersedAdvectionDiffusionWallData3D ************************************ */

template <typename T>
class InstantiateImmersedAdvectionDiffusionWallData3D : public BoxProcessingFunctional3D {
public:
    InstantiateImmersedAdvectionDiffusionWallData3D(
        std::vector<Array<T, 3> > const &vertices_, std::vector<T> const &areas_,
        std::vector<Array<T, 3> > const &normals_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InstantiateImmersedAdvectionDiffusionWallData3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    std::vector<Array<T, 3> > const &vertices;
    std::vector<T> const &areas;
    std::vector<Array<T, 3> > const &normals;
};

template <typename T>
void instantiateImmersedAdvectionDiffusionWallData(
    std::vector<Array<T, 3> > const &vertices, std::vector<T> const &areas,
    MultiContainerBlock3D &container)
{
    static std::vector<Array<T, 3> > dummyNormals;
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    applyProcessingFunctional(
        new InstantiateImmersedAdvectionDiffusionWallData3D<T>(vertices, areas, dummyNormals),
        container.getBoundingBox(), args);
}

template <typename T>
void instantiateImmersedAdvectionDiffusionWallData(
    std::vector<Array<T, 3> > const &vertices, std::vector<T> const &areas,
    std::vector<Array<T, 3> > const &normals, MultiContainerBlock3D &container)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    applyProcessingFunctional(
        new InstantiateImmersedAdvectionDiffusionWallData3D<T>(vertices, areas, normals),
        container.getBoundingBox(), args);
}

/* ******** InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D
 * ************************************ */

template <typename T>
class InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D :
    public BoxProcessingFunctional3D {
public:
    InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D(
        std::vector<Array<T, 3> > const &vertices_, std::vector<T> const &areas_,
        int solutionDomainFlag_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    std::vector<Array<T, 3> > const &vertices;
    std::vector<T> const &areas;
    int solutionDomainFlag;
};

template <typename T>
void instantiateImmersedAdvectionDiffusionWallDataWithTagging(
    std::vector<Array<T, 3> > const &vertices, std::vector<T> const &areas,
    MultiContainerBlock3D &container, MultiScalarField3D<int> &flags, int solutionDomainFlag)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    args.push_back(&flags);
    applyProcessingFunctional(
        new InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D<T>(
            vertices, areas, solutionDomainFlag),
        container.getBoundingBox(), args);
}

/* ******** InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D
 * ************************************ */

// This is the same as InstantiateImmersedAdvectionDiffusionWallDataWithTagging3D, but instead of a
// MultiScalarField3D of flags given to compute the flags of the ImmersedWallData3D, a vector of
// flags on every vertex of the immersed walls is provided.
template <typename T>
class InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D :
    public BoxProcessingFunctional3D {
public:
    InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D(
        std::vector<Array<T, 3> > const &vertices_, std::vector<T> const &areas_,
        std::vector<int> const &flags_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    std::vector<Array<T, 3> > const &vertices;
    std::vector<T> const &areas;
    std::vector<int> const &flags;
};

template <typename T>
void instantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging(
    std::vector<Array<T, 3> > const &vertices, std::vector<T> const &areas,
    std::vector<int> const &flags, MultiContainerBlock3D &container)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    applyProcessingFunctional(
        new InstantiateImmersedAdvectionDiffusionWallDataWithIndexedTagging3D<T>(
            vertices, areas, flags),
        container.getBoundingBox(), args);
}

/* ******** BoxAdvectionDiffusionRhoBarJfunctional3D ************************************ */

template <
    typename T, template <typename U> class Descriptor, template <typename V> class AD_Descriptor>
class BoxAdvectionDiffusionRhoBarJfunctional3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual BoxAdvectionDiffusionRhoBarJfunctional3D<T, Descriptor, AD_Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

}  // namespace plb

#endif  // IMMERSED_ADVECTION_DIFFUSION_WALLS_3D_H
