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

#ifndef OFF_LATTICE_MODEL_3D_H
#define OFF_LATTICE_MODEL_3D_H

#include "atomicBlock/atomicContainerBlock3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"
#include "offLattice/boundaryShapes3D.h"

namespace plb {

template <typename T, class SurfaceData>
class OffLatticeModel3D {
public:
    OffLatticeModel3D(BoundaryShape3D<T, SurfaceData> *shape_, int flowType_);
    OffLatticeModel3D(OffLatticeModel3D<T, SurfaceData> const &rhs);
    OffLatticeModel3D<T, SurfaceData> &operator=(OffLatticeModel3D<T, SurfaceData> const &rhs);
    virtual ~OffLatticeModel3D();
    int getFlowType() const
    {
        return flowType;
    }
    void provideShapeArguments(std::vector<AtomicBlock3D *> args);
    plint getTag(plint id) const;
    bool pointOnSurface(
        Dot3D const &fromPoint, Dot3D const &direction, Array<T, 3> &locatedPoint, T &distance,
        Array<T, 3> &wallNormal, SurfaceData &surfaceData, OffBoundary::Type &bdType,
        plint &id) const;
    Array<T, 3> computeContinuousNormal(
        Array<T, 3> const &p, plint id, bool isAreaWeighted = false) const;
    bool intersectsSurface(Dot3D const &p1, Dot3D const &p2, plint &id) const;
    bool isFluid(Dot3D const &location) const;
    bool isSolid(Dot3D const &location) const;

    bool velIsJ() const
    {
        return velIsJflag;
    }
    void setVelIsJ(bool velIsJflag_)
    {
        velIsJflag = velIsJflag_;
    }
    bool getPartialReplace() const
    {
        return partialReplaceFlag;
    }
    void setPartialReplace(bool prFlag)
    {
        partialReplaceFlag = prFlag;
    }
    void selectSecondOrder(bool flag)
    {
        secondOrderFlag = flag;
    }
    bool usesSecondOrder() const
    {
        return secondOrderFlag;
    }
    void selectUseRegularizedModel(bool flag)
    {
        regularizedModel = flag;
    }
    bool usesRegularizedModel() const
    {
        return regularizedModel;
    }
    void selectComputeStat(bool flag)
    {
        computeStat = flag;
    }
    bool computesStat() const
    {
        return computeStat;
    }
    void setDefineVelocity(bool defineVelocity_)
    {
        defineVelocity = defineVelocity_;
    }
    bool getDefineVelocity() const
    {
        return defineVelocity;
    }

    virtual OffLatticeModel3D<T, SurfaceData> *clone() const = 0;
    virtual plint getNumNeighbors() const = 0;
    virtual bool isExtrapolated() const = 0;
    virtual void prepareCell(Dot3D const &cellLocation, AtomicContainerBlock3D &container) = 0;
    virtual void boundaryCompletion(
        AtomicBlock3D &lattice, AtomicContainerBlock3D &container,
        std::vector<AtomicBlock3D *> const &args) = 0;
    virtual ContainerBlockData *generateOffLatticeInfo() const = 0;
    virtual Array<T, 3> getLocalForce(AtomicContainerBlock3D &container) const = 0;

private:
    BoundaryShape3D<T, SurfaceData> *shape;
    int flowType;
    bool velIsJflag;
    bool partialReplaceFlag;
    bool secondOrderFlag;
    bool regularizedModel;
    bool computeStat;
    bool defineVelocity;
};

/// Precompute a list of nodes which are close to the off-lattice boundary
///   (both wet and dry ones).
template <typename T, class SurfaceData>
class OffLatticePatternFunctional3D : public BoxProcessingFunctional3D {
public:
    /// numNeighbors is the number of neighbors a boundary node must be able
    ///   to access along a given direction.
    OffLatticePatternFunctional3D(OffLatticeModel3D<T, SurfaceData> *offLatticeModel_);
    virtual ~OffLatticePatternFunctional3D();
    OffLatticePatternFunctional3D(OffLatticePatternFunctional3D const &rhs);
    OffLatticePatternFunctional3D &operator=(OffLatticePatternFunctional3D const &rhs);
    void swap(OffLatticePatternFunctional3D &rhs);
    virtual OffLatticePatternFunctional3D<T, SurfaceData> *clone() const;

    /// First AtomicBlock: OffLatticeInfo.
    ///   If there are more atomic-blocks then they are forwarded to the
    ///   shape function, to provide additional read-only parameters.
    /// It is very important that the "offLatticePattern" container block
    /// (the first atomic-block passed) has the same multi-block management
    /// as the lattice used in the simulation.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    OffLatticeModel3D<T, SurfaceData> *offLatticeModel;
};

template <typename T, template <typename U> class Descriptor, class SurfaceData>
class OffLatticeCompletionFunctional3D : public BoxProcessingFunctional3D {
public:
    OffLatticeCompletionFunctional3D(
        OffLatticeModel3D<T, SurfaceData> *offLatticeModel_, plint numShapeArgs_,
        plint numCompletionArgs_);
    virtual ~OffLatticeCompletionFunctional3D();
    OffLatticeCompletionFunctional3D(
        OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData> const &rhs);
    OffLatticeCompletionFunctional3D &operator=(
        OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData> const &rhs);
    void swap(OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData> &rhs);

    /// First AtomicBlock: Lattice; Second AtomicBlock: Container for off-
    ///   lattice info.
    ///   If there are more atomic-blocks then they are forwarded to the
    ///   shape function, to provide additional read-only parameters.
    ///
    /// It is very important that the "offLatticePattern" container block
    /// which is passed as the second atomic-block with the off-lattice info
    /// has the same multi-block management as the lattice used in the simulation.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual OffLatticeCompletionFunctional3D<T, Descriptor, SurfaceData> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    OffLatticeModel3D<T, SurfaceData> *offLatticeModel;
    plint numShapeArgs, numCompletionArgs;
};

template <typename T, class SurfaceData>
class GetForceOnObjectFunctional3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    GetForceOnObjectFunctional3D(OffLatticeModel3D<T, SurfaceData> *offLatticeModel_);
    virtual ~GetForceOnObjectFunctional3D();
    GetForceOnObjectFunctional3D(GetForceOnObjectFunctional3D<T, SurfaceData> const &rhs);
    GetForceOnObjectFunctional3D<T, SurfaceData> &operator=(
        GetForceOnObjectFunctional3D<T, SurfaceData> const &rhs);

    /// First AtomicBlock: Container for off-lattice info.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual GetForceOnObjectFunctional3D<T, SurfaceData> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    Array<T, 3> getForce() const;

private:
    OffLatticeModel3D<T, SurfaceData> *offLatticeModel;
    Array<plint, 3> forceId;
};

template <typename T, class BoundaryType>
Array<T, 3> getForceOnObject(
    MultiBlock3D &offLatticePattern, OffLatticeModel3D<T, BoundaryType> const &offLatticeModel);

}  // namespace plb

#endif  // OFF_LATTICE_MODEL_3D_H
