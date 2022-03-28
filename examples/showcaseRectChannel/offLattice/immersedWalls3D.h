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

#ifndef IMMERSED_WALLS_3D_H
#define IMMERSED_WALLS_3D_H

#include <memory>

#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "multiBlock/multiBlockGenerator3D.h"

namespace plb {

/* ******** FavierDeltaFunction ************************************ */

template <typename T>
class FavierDeltaFunction {
public:
    FavierDeltaFunction(int N_) : N(N_), dx(4. / (T)N), invDx(1. / dx)
    {
        sampleFunction();
    }
    T rawValue(T r) const
    {
        T rabs = std::fabs(r);
        T rsqr = r * r;
        if (rabs < 0.5) {
            return (1 + sqrt(-3.0 * rsqr + 1.0)) / 3.0;
        } else if (rabs < 1.5) {
            return (5.0 - 3.0 * rabs - sqrt(-3.0 * (1 - rabs) * (1.0 - rabs) + 1.0)) / 6.0;
        } else {
            return 0.;
        }
    }
    T w(T r) const
    {
        int position = (int)((r + 2.0) * invDx + 0.5);
        if (position <= 0) {
            return 0.;
        }
        if (position >= N) {
            return 0.;
        }
        return samples[position];
    }
    T W(Array<T, 3> const &r) const
    {
        return w(r[0]) * w(r[1]) * w(r[2]);
    }

private:
    void sampleFunction()
    {
        samples.resize(N + 1);
        for (int i = 0; i <= N; ++i) {
            samples[i] = rawValue(-2. + dx * i);
        }
    }

private:
    int N;
    T dx, invDx;
    std::vector<T> samples;
};

template <typename T>
FavierDeltaFunction<T> const &favierDeltaFunction()
{
    static FavierDeltaFunction<T> deltaFunction(1000);
    return deltaFunction;
}

/* ******** InamuroDeltaFunction ************************************ */

template <typename T>
class InamuroDeltaFunction {
public:
    InamuroDeltaFunction(int N_) : N(N_), dx(4. / (T)N), invDx(1. / dx)
    {
        sampleFunction();
    }
    T rawValue(T r) const
    {
        T rabs = std::fabs(r);
        T rsqr = r * r;
        if (rabs < 1.) {
            return 0.125 * (3. - 2. * rabs + std::sqrt(1. + 4. * rabs - 4. * rsqr));
        } else if (rabs < 2.) {
            return 0.125 * (5. - 2. * rabs - std::sqrt(-7. + 12. * rabs - 4. * rsqr));
        } else {
            return 0.;
        }
    }
    T w(T r) const
    {
        int position = (int)((r + 2.0) * invDx + 0.5);
        if (position <= 0) {
            return 0.;
        }
        if (position >= N) {
            return 0.;
        }
        return samples[position];
    }
    T W(Array<T, 3> const &r) const
    {
        return w(r[0]) * w(r[1]) * w(r[2]);
    }

private:
    void sampleFunction()
    {
        samples.resize(N + 1);
        for (int i = 0; i <= N; ++i) {
            samples[i] = rawValue(-2. + dx * i);
        }
    }

private:
    int N;
    T dx, invDx;
    std::vector<T> samples;
};

template <typename T>
inline InamuroDeltaFunction<T> const &inamuroDeltaFunction()
{
    static InamuroDeltaFunction<T> deltaFunction(1000);
    return deltaFunction;
}

/* ******** ImmersedWallData3D ************************************ */

template <typename T>
struct ImmersedWallData3D : public ContainerBlockData {
    std::vector<Array<T, 3> > vertices;
    std::vector<T> areas;
    std::vector<Array<T, 3> > normals;
    std::vector<Array<T, 3> > g;
    std::vector<int> flags;  // Flag for each vertex used to distinguish between vertices for
                             // conditional reduction operations.
    std::vector<pluint> globalVertexIds;
    virtual ImmersedWallData3D<T> *clone() const
    {
        return new ImmersedWallData3D<T>(*this);
    }
};

/* ******** SurfaceBlockData3D ************************************ */

template <typename T>
struct SurfaceBlockData3D : public ContainerBlockData {
    std::vector<Array<T, 3> > vertices;

    virtual SurfaceBlockData3D<T> *clone() const
    {
        return new SurfaceBlockData3D<T>(*this);
    }
};

/* ******** Utility functions ************************************ */

template <typename T>
inline bool closedOpenContained(Array<T, 3> const &x, Box3D const &box)
{
    return x[0] >= (box.x0 - 0.5) && x[0] < (box.x1 + 0.5) && x[1] >= (box.y0 - 0.5)
           && x[1] < (box.y1 + 0.5) && x[2] >= (box.z0 - 0.5) && x[2] < (box.z1 + 0.5);
    // in order to count correctly the particles, a 0.5 must be added
}

/* ******** ReduceImmersedTorque3D ************************************ */

// The reduced quantity is computed only for the vertices which have a flag
// equal to "reductionFlag".
template <typename T>
class ReduceImmersedTorque3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    ReduceImmersedTorque3D(Array<T, 3> const &center_, int reductionFlag_ = 0);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ReduceImmersedTorque3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    Array<T, 3> getSumTorque() const;

private:
    Array<T, 3> center;
    Array<plint, 3> sum_torque_ids;
    int reductionFlag;
};

template <typename T>
Array<T, 3> reduceImmersedTorque(
    MultiContainerBlock3D &container, Array<T, 3> const &center, int reductionFlag = 0)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    ReduceImmersedTorque3D<T> functional(center, reductionFlag);
    applyProcessingFunctional(functional, container.getBoundingBox(), args);
    return functional.getSumTorque();
}

/* ******** ReduceImmersedForce3D ************************************ */

// The reduced quantity is computed only for the vertices which have a flag
// equal to "reductionFlag".
template <typename T>
class ReduceImmersedForce3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    ReduceImmersedForce3D(int reductionFlag_ = 0);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ReduceImmersedForce3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    Array<T, 3> getSumG() const;

private:
    Array<plint, 3> sum_g_ids;
    int reductionFlag;
};

template <typename T>
Array<T, 3> reduceImmersedForce(MultiContainerBlock3D &container, int reductionFlag = 0)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    ReduceImmersedForce3D<T> functional(reductionFlag);
    applyProcessingFunctional(functional, container.getBoundingBox(), args);
    return functional.getSumG();
}

/* ******** ReduceImmersedArea3D ************************************ */

// The reduced quantity is computed only for the vertices which have a flag
// equal to "reductionFlag".
template <typename T>
class ReduceImmersedArea3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    ReduceImmersedArea3D(int reductionFlag_ = 0);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ReduceImmersedArea3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    T getSumArea() const;

private:
    plint sum_area_id;
    int reductionFlag;
};

template <typename T>
T reduceImmersedArea(MultiContainerBlock3D &container, int reductionFlag = 0)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    ReduceImmersedArea3D<T> functional(reductionFlag);
    applyProcessingFunctional(functional, container.getBoundingBox(), args);
    return functional.getSumArea();
}

/* ******** InamuroIteration3D ************************************ */

template <typename T, class VelFunction>
class InamuroIteration3D : public BoxProcessingFunctional3D {
public:
    InamuroIteration3D(VelFunction velFunction_, T tau_, bool incompressibleModel_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InamuroIteration3D<T, VelFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    VelFunction velFunction;
    T tau;
    bool incompressibleModel;
};

template <typename T, class VelFunction>
void inamuroIteration(
    VelFunction velFunction, MultiScalarField3D<T> &rhoBar, MultiTensorField3D<T, 3> &j,
    MultiContainerBlock3D &container, T tau, bool incompressibleModel)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&rhoBar);
    args.push_back(&j);
    args.push_back(&container);
    applyProcessingFunctional(
        new InamuroIteration3D<T, VelFunction>(velFunction, tau, incompressibleModel),
        rhoBar.getBoundingBox(), args);
}

/* ******** IndexedInamuroIteration3D ************************************ */

// This is the same as InamuroIteration3D, with the difference that
// the VelFunction accepts as argument a global vertex index instead of
// a 3D position in space.
template <typename T, class VelFunction>
class IndexedInamuroIteration3D : public BoxProcessingFunctional3D {
public:
    IndexedInamuroIteration3D(VelFunction velFunction_, T tau_, bool incompressibleModel_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual IndexedInamuroIteration3D<T, VelFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    VelFunction velFunction;
    T tau;
    bool incompressibleModel;
};

template <typename T, class VelFunction>
void indexedInamuroIteration(
    VelFunction velFunction, MultiScalarField3D<T> &rhoBar, MultiTensorField3D<T, 3> &j,
    MultiContainerBlock3D &container, T tau, bool incompressibleModel)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&rhoBar);
    args.push_back(&j);
    args.push_back(&container);
    applyProcessingFunctional(
        new IndexedInamuroIteration3D<T, VelFunction>(velFunction, tau, incompressibleModel),
        rhoBar.getBoundingBox(), args);
}

/* ******** IndexedImmersedBoundaryIteration3D ************************************ */

template <typename T, class VelFunction>
class IndexedImmersedBoundaryIteration3D : public BoxProcessingFunctional3D {
public:
    IndexedImmersedBoundaryIteration3D(VelFunction velFunction_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual IndexedImmersedBoundaryIteration3D<T, VelFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    VelFunction velFunction;
};

template <typename T, class VelFunction>
void indexedImmersedBoundaryIteration(
    VelFunction velFunction, MultiTensorField3D<T, 3> &u, MultiContainerBlock3D &container)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&u);
    args.push_back(&container);
    applyProcessingFunctional(
        new IndexedImmersedBoundaryIteration3D<T, VelFunction>(velFunction), u.getBoundingBox(),
        args);
}

/* ******** ConstVelInamuroIteration3D ************************************ */

template <typename T>
class ConstVelInamuroIteration3D : public BoxProcessingFunctional3D {
public:
    ConstVelInamuroIteration3D(Array<T, 3> const &wallVelocity_, T tau_, bool incompressibleModel_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ConstVelInamuroIteration3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    Array<T, 3> wallVelocity;
    T tau;
    bool incompressibleModel;
};

template <typename T>
void constVelInamuroIteration(
    Array<T, 3> const &wallVelocity, MultiScalarField3D<T> &rhoBar, MultiTensorField3D<T, 3> &j,
    MultiContainerBlock3D &container, T tau, bool incompressibleModel)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&rhoBar);
    args.push_back(&j);
    args.push_back(&container);
    applyProcessingFunctional(
        new ConstVelInamuroIteration3D<T>(wallVelocity, tau, incompressibleModel),
        rhoBar.getBoundingBox(), args);
}

/* ******** ComputeImmersedBoundaryForce3D ************************************ */

// This functional computes the immersed boundary force on the lattice and
// stores it in a provided tensor field. This data processor must be called
// after all the immersed boundary iterations have completed.
template <typename T>
class ComputeImmersedBoundaryForce3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ComputeImmersedBoundaryForce3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
void computeImmersedBoundaryForce3D(
    MultiTensorField3D<T, 3> &force, MultiContainerBlock3D &container)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&force);
    args.push_back(&container);
    applyProcessingFunctional(new ComputeImmersedBoundaryForce3D<T>, force.getBoundingBox(), args);
}

/* ******** InstantiateImmersedWallData3D ************************************ */

template <typename T>
class InstantiateImmersedWallData3D : public BoxProcessingFunctional3D {
public:
    InstantiateImmersedWallData3D(
        std::vector<Array<T, 3> > const &vertices_, std::vector<T> const &areas_,
        std::vector<Array<T, 3> > const &normals_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InstantiateImmersedWallData3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    std::vector<Array<T, 3> > const &vertices;
    std::vector<T> const &areas;
    std::vector<Array<T, 3> > const &normals;
};

template <typename T>
void instantiateImmersedWallData(
    std::vector<Array<T, 3> > const &vertices, std::vector<T> const &areas,
    MultiContainerBlock3D &container)
{
    static std::vector<Array<T, 3> > dummyNormals;
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    applyProcessingFunctional(
        new InstantiateImmersedWallData3D<T>(vertices, areas, dummyNormals),
        container.getBoundingBox(), args);
}

template <typename T>
void instantiateImmersedWallData(
    std::vector<Array<T, 3> > const &vertices, std::vector<T> const &areas,
    std::vector<Array<T, 3> > const &normals, MultiContainerBlock3D &container)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    applyProcessingFunctional(
        new InstantiateImmersedWallData3D<T>(vertices, areas, normals), container.getBoundingBox(),
        args);
}

/* ******** InstantiateImmersedWallDataWithTagging3D ************************************ */

template <typename T>
class InstantiateImmersedWallDataWithTagging3D : public BoxProcessingFunctional3D {
public:
    InstantiateImmersedWallDataWithTagging3D(
        std::vector<Array<T, 3> > const &vertices_, std::vector<T> const &areas_, int fluidFlag_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InstantiateImmersedWallDataWithTagging3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    std::vector<Array<T, 3> > const &vertices;
    std::vector<T> const &areas;
    int fluidFlag;
};

template <typename T>
void instantiateImmersedWallDataWithTagging(
    std::vector<Array<T, 3> > const &vertices, std::vector<T> const &areas,
    MultiContainerBlock3D &container, MultiScalarField3D<int> &flags, int fluidFlag)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    args.push_back(&flags);
    applyProcessingFunctional(
        new InstantiateImmersedWallDataWithTagging3D<T>(vertices, areas, fluidFlag),
        container.getBoundingBox(), args);
}

/* ******** InstantiateImmersedWallDataWithIndexedTagging3D ************************************ */

// This is the same as InstantiateImmersedWallDataWithTagging3D, but instead of a MultiScalarField3D
// of flags given to compute the flags of the ImmersedWallData3D, a vector of flags on every
// vertex of the immersed walls is provided.
template <typename T>
class InstantiateImmersedWallDataWithIndexedTagging3D : public BoxProcessingFunctional3D {
public:
    InstantiateImmersedWallDataWithIndexedTagging3D(
        std::vector<Array<T, 3> > const &vertices_, std::vector<T> const &areas_,
        std::vector<int> const &flags_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InstantiateImmersedWallDataWithIndexedTagging3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    std::vector<Array<T, 3> > const &vertices;
    std::vector<T> const &areas;
    std::vector<int> const &flags;
};

template <typename T>
void instantiateImmersedWallDataWithIndexedTagging(
    std::vector<Array<T, 3> > const &vertices, std::vector<T> const &areas,
    std::vector<int> const &flags, MultiContainerBlock3D &container)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    applyProcessingFunctional(
        new InstantiateImmersedWallDataWithIndexedTagging3D<T>(vertices, areas, flags),
        container.getBoundingBox(), args);
}

/* ******** InstantiateSurfaceBlockData3D ****************************** */

template <typename T>
class InstantiateSurfaceBlockData3D : public BoxProcessingFunctional3D {
public:
    InstantiateSurfaceBlockData3D(plint envelopeWidth_, std::vector<Array<T, 3> > const &vertices_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual InstantiateSurfaceBlockData3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    plint envelopeWidth;
    std::vector<Array<T, 3> > const &vertices;
};

template <typename T>
void instantiateSurfaceBlockData(
    plint envelopeWidth, std::vector<Array<T, 3> > const &vertices,
    MultiContainerBlock3D &container)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    applyProcessingFunctional(
        new InstantiateSurfaceBlockData3D<T>(envelopeWidth, vertices), container.getBoundingBox(),
        args);
}

/* ******** SurfaceOnLattice3D **************************************** */

template <typename T, typename U>
class SurfaceOnLattice3D : public BoxProcessingFunctional3D {
public:
    SurfaceOnLattice3D(U value_, plint envelopeWidth_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual SurfaceOnLattice3D<T, U> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    U value;
    plint envelopeWidth;
};

template <typename T, typename U>
void surfaceOnLattice(
    U value, plint envelopeWidth, MultiScalarField3D<U> &values, MultiContainerBlock3D &container)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&values);
    args.push_back(&container);
    applyProcessingFunctional(
        new SurfaceOnLattice3D<T, U>(value, envelopeWidth), container.getBoundingBox(), args);
}

/* ******** SurfaceOnLattice3D_N **************************************** */

// This is the same as SurfaceOnLattice3D, but instead of the values being
// stored in a MultiScalarField3D<U>, they are stored in a MultiNTensorField3D<U>.

template <typename T, typename U>
class SurfaceOnLattice3D_N : public BoxProcessingFunctional3D {
public:
    SurfaceOnLattice3D_N(U value_, plint envelopeWidth_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual SurfaceOnLattice3D_N<T, U> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    U value;
    plint envelopeWidth;
};

template <typename T, typename U>
void surfaceOnLattice_N(
    U value, plint envelopeWidth, MultiNTensorField3D<U> &values, MultiContainerBlock3D &container)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&values);
    args.push_back(&container);
    applyProcessingFunctional(
        new SurfaceOnLattice3D_N<T, U>(value, envelopeWidth), container.getBoundingBox(), args);
}

/* ******** ResetForceStatistics3D ************************************ */

// This data processor resets to zero the "per surface vertex" force vectors
// which reside in the immersed data container field. This is used for
// optimization purposes. Sometimes when the surface is not moving, the user
// should instantiate the immersed wall data only once, and not at every
// itaration. Doing so, would not work for the force computations, since the
// forces are added up during the Inamuro iterations. This is why, before
// measuring the forces, one must call this data processor, so that the force
// variable in the container is set back to zero.
template <typename T>
class ResetForceStatistics3D : public BoxProcessingFunctional3D {
public:
    ResetForceStatistics3D() { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ResetForceStatistics3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
void resetForceStatistics(MultiContainerBlock3D &container)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&container);
    applyProcessingFunctional(new ResetForceStatistics3D<T>(), container.getBoundingBox(), args);
}

/* ******** RecomputeImmersedForce3D ************************************ */

// This class recomputes the immersed force (variable "g" in the ImmersedWallData3D)
// by using the classical stress tensor relation.
// The normalFunction is a function:
//   Array<T,3> normalFunction(plint id);
// which takes a global vertex id and computes the unit normal at that point.
template <typename T, template <typename U> class Descriptor, class NormalFunction>
class RecomputeImmersedForce3D : public BoxProcessingFunctional3D {
public:
    RecomputeImmersedForce3D(
        NormalFunction normalFunction_, T omega_, T densityOffset_, bool incompressibleModel_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual RecomputeImmersedForce3D<T, Descriptor, NormalFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    NormalFunction normalFunction;
    T omega;
    T rho0;
    bool incompressibleModel;
};

template <typename T, template <typename U> class Descriptor, class NormalFunction>
void recomputeImmersedForce(
    NormalFunction normalFunction, T omega, T densityOffset,
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiContainerBlock3D &container,
    plint envelopeWidth, Box3D domain, bool incompressibleModel)
{
    std::unique_ptr<MultiScalarField3D<T> > rhoBar =
        generateMultiScalarField<T>(lattice, envelopeWidth);
    computeRhoBar<T, Descriptor>(lattice, *rhoBar, domain);

    std::unique_ptr<MultiTensorField3D<T, SymmetricTensorImpl<T, 3>::n> > PiNeq =
        generateMultiTensorField<T, SymmetricTensorImpl<T, 3>::n>(lattice, envelopeWidth);
    computePiNeq<T, Descriptor>(lattice, *PiNeq, domain);

    std::vector<MultiBlock3D *> args;
    args.push_back(rhoBar.get());
    args.push_back(PiNeq.get());
    args.push_back(&container);
    applyProcessingFunctional(
        new RecomputeImmersedForce3D<T, Descriptor, NormalFunction>(
            normalFunction, omega, densityOffset, incompressibleModel),
        domain, args);
}

/* ******** OpenSurfaceImmersedForce3D ************************************ */

// This class recomputes the immersed force (variable "g" in the ImmersedWallData3D)
// by using the classical stress tensor relation.
// The normalFunction is a function:
//   Array<T,3> normalFunction(plint id);
// which takes a global vertex id and computes the unit normal at that point.
template <typename T, template <typename U> class Descriptor, class NormalFunction>
class OpenSurfaceImmersedForce3D : public BoxProcessingFunctional3D {
public:
    OpenSurfaceImmersedForce3D(
        NormalFunction normalFunction_, T omega_, T densityOffset_, bool incompressibleModel_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual OpenSurfaceImmersedForce3D<T, Descriptor, NormalFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    NormalFunction normalFunction;
    T omega;
    T rho0;
    bool incompressibleModel;
};

template <typename T, template <typename U> class Descriptor, class NormalFunction>
void openSurfaceImmersedForce(
    NormalFunction normalFunction, T omega, T densityOffset,
    MultiBlockLattice3D<T, Descriptor> &lattice, MultiContainerBlock3D &container,
    plint envelopeWidth, Box3D domain, bool incompressibleModel)
{
    std::unique_ptr<MultiScalarField3D<T> > rhoBar =
        generateMultiScalarField<T>(lattice, envelopeWidth);
    computeRhoBar<T, Descriptor>(lattice, *rhoBar, domain);

    std::unique_ptr<MultiTensorField3D<T, SymmetricTensorImpl<T, 3>::n> > PiNeq =
        generateMultiTensorField<T, SymmetricTensorImpl<T, 3>::n>(lattice, envelopeWidth);
    computePiNeq<T, Descriptor>(lattice, *PiNeq, domain);

    std::vector<MultiBlock3D *> args;
    args.push_back(rhoBar.get());
    args.push_back(PiNeq.get());
    args.push_back(&container);
    applyProcessingFunctional(
        new OpenSurfaceImmersedForce3D<T, Descriptor, NormalFunction>(
            normalFunction, omega, densityOffset, incompressibleModel),
        domain, args);
}

/* ******** TwoPhaseInamuroParam3D ************************************ */

template <typename T>
class TwoPhaseInamuroParam3D {
public:
    TwoPhaseInamuroParam3D(std::vector<AtomicBlock3D *> &blocks, T tau_, T tau2_);
    pluint getNumVertices() const
    {
        return numVertices;
    }
    T area(plint i) const;
    Array<T, 3> &g(plint i);
    Array<T, 3> vertex(plint i) const;
    T rhoBar(plint iX, plint iY, plint iZ) const;
    Array<T, 3> j(plint iX, plint iY, plint iZ) const;
    void addToJ(plint iX, plint iY, plint iZ, Array<T, 3> deltaJ);
    T getTau(plint iX, plint iY, plint iZ) const;
    Dot3D getLocation() const;

private:
    int getFlag(plint iX, plint iY, plint iZ) const;
    pluint getGlobalVertexId(plint i) const;

private:
    ScalarField3D<T> *rhoBar_, *rhoBar2_;
    TensorField3D<T, 3> *j_, *j2_;
    ScalarField3D<int> *flag_;
    ScalarField3D<T> *volumeFraction_;
    AtomicContainerBlock3D *container;
    ImmersedWallData3D<T> *wallData;
    Dot3D ofsRhoBar2, ofsJ, ofsJ2, ofsFlag, ofsVF;
    Dot3D location;
    pluint numVertices;
    T tau, tau2;
};

/* ******** TwoPhaseInamuroIteration3D ************************************ */

template <typename T, class VelFunction>
class TwoPhaseInamuroIteration3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseInamuroIteration3D(VelFunction velFunction_, T tau_, T tau2_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual TwoPhaseInamuroIteration3D<T, VelFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    VelFunction velFunction;
    T tau, tau2;
};

template <typename T, class VelFunction>
void twoPhaseInamuroIteration(
    VelFunction velFunction, MultiScalarField3D<T> &rhoBar, MultiScalarField3D<T> &rhoBar2,
    MultiTensorField3D<T, 3> &j, MultiTensorField3D<T, 3> &j2, MultiScalarField3D<int> &flag,
    MultiScalarField3D<T> &volumeFraction, MultiContainerBlock3D &container, T tau, T tau2)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&rhoBar);
    args.push_back(&rhoBar2);
    args.push_back(&j);
    args.push_back(&j2);
    args.push_back(&flag);
    args.push_back(&volumeFraction);
    args.push_back(&container);
    applyProcessingFunctional(
        new TwoPhaseInamuroIteration3D<T, VelFunction>(velFunction, tau, tau2),
        rhoBar.getBoundingBox(), args);
}

/* ******** TwoPhaseIndexedInamuroIteration3D ************************************ */

// This is the same as TwoPhaseInamuroIteration3D, with the difference that
// the VelFunction accepts as argument a global vertex index instead of
// a 3D position in space.
template <typename T, class VelFunction>
class TwoPhaseIndexedInamuroIteration3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseIndexedInamuroIteration3D(VelFunction velFunction_, T tau_, T tau2_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual TwoPhaseIndexedInamuroIteration3D<T, VelFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    VelFunction velFunction;
    T tau, tau2;
};

template <typename T, class VelFunction>
void twoPhaseIndexedInamuroIteration(
    VelFunction velFunction, MultiScalarField3D<T> &rhoBar, MultiScalarField3D<T> &rhoBar2,
    MultiTensorField3D<T, 3> &j, MultiTensorField3D<T, 3> &j2, MultiScalarField3D<int> &flag,
    MultiScalarField3D<T> &volumeFraction, MultiContainerBlock3D &container, T tau, T tau2)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&rhoBar);
    args.push_back(&rhoBar2);
    args.push_back(&j);
    args.push_back(&j2);
    args.push_back(&flag);
    args.push_back(&volumeFraction);
    args.push_back(&container);
    applyProcessingFunctional(
        new TwoPhaseIndexedInamuroIteration3D<T, VelFunction>(velFunction, tau, tau2),
        rhoBar.getBoundingBox(), args);
}

/* ******** TwoPhaseConstVelInamuroIteration3D ************************************ */

template <typename T>
class TwoPhaseConstVelInamuroIteration3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseConstVelInamuroIteration3D(Array<T, 3> const &wallVelocity_, T tau_, T tau2_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual TwoPhaseConstVelInamuroIteration3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    Array<T, 3> wallVelocity;
    T tau, tau2;
};

template <typename T>
void twoPhaseConstVelInamuroIteration(
    Array<T, 3> const &wallVelocity, MultiScalarField3D<T> &rhoBar, MultiScalarField3D<T> &rhoBar2,
    MultiTensorField3D<T, 3> &j, MultiTensorField3D<T, 3> &j2, MultiScalarField3D<int> &flag,
    MultiScalarField3D<T> &volumeFraction, MultiContainerBlock3D &container, T tau, T tau2)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&rhoBar);
    args.push_back(&rhoBar2);
    args.push_back(&j);
    args.push_back(&j2);
    args.push_back(&flag);
    args.push_back(&volumeFraction);
    args.push_back(&container);
    applyProcessingFunctional(
        new TwoPhaseConstVelInamuroIteration3D<T>(wallVelocity, tau, tau2), rhoBar.getBoundingBox(),
        args);
}

}  // namespace plb

#endif  // IMMERSED_WALLS_3D_H
