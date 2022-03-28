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

#ifndef GENERALIZED_OFF_LATTICE_MODEL_3D_H
#define GENERALIZED_OFF_LATTICE_MODEL_3D_H

#include "core/globalDefs.h"
#include "offLattice/offLatticeModel3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class ExtrapolatedGeneralizedOffLatticeModel3D : public OffLatticeModel3D<T, Array<T, 3> > {
public:
    ExtrapolatedGeneralizedOffLatticeModel3D(
        BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_);
    virtual ExtrapolatedGeneralizedOffLatticeModel3D<T, Descriptor> *clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const;
    virtual void prepareCell(Dot3D const &cellLocation, AtomicContainerBlock3D &container);
    virtual void boundaryCompletion(
        AtomicBlock3D &lattice, AtomicContainerBlock3D &container,
        std::vector<AtomicBlock3D *> const &args);
    virtual ContainerBlockData *generateOffLatticeInfo() const;
    virtual Array<T, 3> getLocalForce(AtomicContainerBlock3D &container) const;

private:
    void cellCompletion(
        BlockLattice3D<T, Descriptor> &lattice, Dot3D const &genNode,
        std::vector<std::pair<int, int> > const &dryNodeFluidDirections,
        std::vector<int> const &dryNodeFluidNoSolidDirections, std::vector<plint> const &dryNodeIds,
        Dot3D const &absoluteOffset, Array<T, 3> &localForce,
        std::vector<AtomicBlock3D *> const &args);

private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class ExtrapolatedGeneralizedOffLatticeInfo3D : public ContainerBlockData {
    public:
        ExtrapolatedGeneralizedOffLatticeInfo3D() : localForce(Array<T, 3>::zero()) { }
        std::vector<Dot3D> const &getDryNodes() const
        {
            return dryNodes;
        }
        std::vector<Dot3D> &getDryNodes()
        {
            return dryNodes;
        }
        std::vector<std::vector<std::pair<int, int> > > const &getDryNodeFluidDirections() const
        {
            return dryNodeFluidDirections;
        }
        std::vector<std::vector<std::pair<int, int> > > &getDryNodeFluidDirections()
        {
            return dryNodeFluidDirections;
        }
        std::vector<std::vector<int> > const &getDryNodeFluidWithFluidDirections() const
        {
            return dryNodeFluidNoSolidDirections;
        }
        std::vector<std::vector<int> > &getDryNodeFluidWithFluidDirections()
        {
            return dryNodeFluidNoSolidDirections;
        }
        std::vector<std::vector<plint> > const &getDryNodeIds() const
        {
            return dryNodeIds;
        }
        std::vector<std::vector<plint> > &getDryNodeIds()
        {
            return dryNodeIds;
        }
        Array<T, 3> const &getLocalForce() const
        {
            return localForce;
        }
        Array<T, 3> &getLocalForce()
        {
            return localForce;
        }
        virtual ExtrapolatedGeneralizedOffLatticeInfo3D *clone() const
        {
            return new ExtrapolatedGeneralizedOffLatticeInfo3D(*this);
        }

    private:
        std::vector<Dot3D> dryNodes;
        std::vector<std::vector<std::pair<int, int> > > dryNodeFluidDirections;
        std::vector<std::vector<int> > dryNodeFluidNoSolidDirections;
        std::vector<std::vector<plint> > dryNodeIds;
        Array<T, 3> localForce;
    };
};

template <typename T, template <typename U> class Descriptor>
class InterpolatedGeneralizedOffLatticeModel3D : public OffLatticeModel3D<T, Array<T, 3> > {
public:
    InterpolatedGeneralizedOffLatticeModel3D(
        BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_);
    virtual InterpolatedGeneralizedOffLatticeModel3D<T, Descriptor> *clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const;
    virtual void prepareCell(Dot3D const &cellLocation, AtomicContainerBlock3D &container);
    virtual void boundaryCompletion(
        AtomicBlock3D &lattice, AtomicContainerBlock3D &container,
        std::vector<AtomicBlock3D *> const &args);
    virtual ContainerBlockData *generateOffLatticeInfo() const;
    virtual Array<T, 3> getLocalForce(AtomicContainerBlock3D &container) const;

private:
    void cellCompletion(
        BlockLattice3D<T, Descriptor> &lattice, Dot3D const &genNode,
        std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections,
        std::vector<int> const &wetNodeFluidDirections, std::vector<plint> const &wetNodeIds,
        std::vector<plint> const &allWetNodeIds, std::vector<int> const &solidDirections,
        Dot3D const &absoluteOffset, Array<T, 3> &localForce,
        std::vector<AtomicBlock3D *> const &args);

private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class InterpolatedGeneralizedOffLatticeInfo3D : public ContainerBlockData {
    public:
        std::vector<Dot3D> const &getWetNodes() const
        {
            return wetNodes;
        }
        std::vector<Dot3D> &getWetNodes()
        {
            return wetNodes;
        }
        std::vector<std::vector<std::pair<int, int> > > const &getWetNodeSolidUsableDirections()
            const
        {
            return wetNodeSolidUsableDirections;
        }
        std::vector<std::vector<std::pair<int, int> > > &getWetNodeSolidUsableDirections()
        {
            return wetNodeSolidUsableDirections;
        }
        std::vector<std::vector<int> > const &getWetNodeFluidDirections() const
        {
            return wetNodeFluidDirections;
        }
        std::vector<std::vector<int> > &getWetNodeFluidDirections()
        {
            return wetNodeFluidDirections;
        }
        std::vector<std::vector<int> > const &getSolidNeighbors() const
        {
            return solidNeighbors;
        }
        std::vector<std::vector<int> > &getSolidNeighbors()
        {
            return solidNeighbors;
        }
        std::vector<std::vector<plint> > const &getWetNodeIds() const
        {
            return wetNodeIds;
        }
        std::vector<std::vector<plint> > &getWetNodeIds()
        {
            return wetNodeIds;
        }
        std::vector<std::vector<plint> > const &getAllWetNodeIds() const
        {
            return allWetNodeIds;
        }
        std::vector<std::vector<plint> > &getAllWetNodeIds()
        {
            return allWetNodeIds;
        }
        Array<T, 3> const &getLocalForce() const
        {
            return localForce;
        }
        Array<T, 3> &getLocalForce()
        {
            return localForce;
        }
        virtual InterpolatedGeneralizedOffLatticeInfo3D *clone() const
        {
            return new InterpolatedGeneralizedOffLatticeInfo3D(*this);
        }

    private:
        std::vector<Dot3D> wetNodes;
        std::vector<std::vector<std::pair<int, int> > >
            wetNodeSolidUsableDirections;  // stores the directions where there is a wall and valid
                                           // neighbors in its opposite direction
        std::vector<std::vector<int> > wetNodeFluidDirections;  //
        std::vector<std::vector<int> > solidNeighbors;          //
        std::vector<std::vector<plint> > wetNodeIds, allWetNodeIds;
        Array<T, 3> localForce;
    };
};

template <typename T, template <typename U> class Descriptor>
class InterpolatedFdOffLatticeModel3D : public OffLatticeModel3D<T, Array<T, 3> > {
public:
    InterpolatedFdOffLatticeModel3D(BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_);
    virtual InterpolatedFdOffLatticeModel3D<T, Descriptor> *clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const;
    virtual void prepareCell(Dot3D const &cellLocation, AtomicContainerBlock3D &container);
    virtual void boundaryCompletion(
        AtomicBlock3D &lattice, AtomicContainerBlock3D &container,
        std::vector<AtomicBlock3D *> const &args);
    virtual ContainerBlockData *generateOffLatticeInfo() const;
    virtual Array<T, 3> getLocalForce(AtomicContainerBlock3D &container) const;

private:
    void cellCompletion(
        BlockLattice3D<T, Descriptor> &lattice, Dot3D const &genNode,
        std::vector<std::pair<int, int> > const &wetNodeSolidUsableDirections,
        std::vector<int> const &wetNodeFluidDirections, std::vector<plint> const &wetNodeIds,
        std::vector<int> const &solidDirections, Dot3D const &absoluteOffset,
        const std::pair<int, int> &xDerivDirAndOrder, const std::pair<int, int> &yDerivDirAndOrder,
        const std::pair<int, int> &zDerivDirAndOrder, Array<T, 3> &localForce,
        std::vector<AtomicBlock3D *> const &args);
    bool isUsable(const Dot3D &pos) const;
    std::pair<int, int> computeOrderAndDirection(const Dot3D &pos, const Dot3D &dx) const;

private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class InterpolatedFdOffLatticeInfo3D : public ContainerBlockData {
    public:
        std::vector<Dot3D> const &getWetNodes() const
        {
            return wetNodes;
        }
        std::vector<Dot3D> &getWetNodes()
        {
            return wetNodes;
        }
        std::vector<std::vector<std::pair<int, int> > > const &getWetNodeSolidUsableDirections()
            const
        {
            return wetNodeSolidUsableDirections;
        }
        std::vector<std::vector<std::pair<int, int> > > &getWetNodeSolidUsableDirections()
        {
            return wetNodeSolidUsableDirections;
        }
        std::vector<std::vector<int> > const &getWetNodeFluidDirections() const
        {
            return wetNodeFluidDirections;
        }
        std::vector<std::vector<int> > &getWetNodeFluidDirections()
        {
            return wetNodeFluidDirections;
        }
        std::vector<std::vector<plint> > const &getWetNodeIds() const
        {
            return wetNodeIds;
        }
        std::vector<std::vector<plint> > &getWetNodeIds()
        {
            return wetNodeIds;
        }
        std::vector<std::vector<int> > const &getSolidNeighbors() const
        {
            return solidNeighbors;
        }
        std::vector<std::vector<int> > &getSolidNeighbors()
        {
            return solidNeighbors;
        }
        std::vector<std::pair<int, int> > const &getXderivDirAndOrder() const
        {
            return xDerivDirAndOrder;
        }
        std::vector<std::pair<int, int> > &getXderivDirAndOrder()
        {
            return xDerivDirAndOrder;
        }
        std::vector<std::pair<int, int> > const &getYderivDirAndOrder() const
        {
            return yDerivDirAndOrder;
        }
        std::vector<std::pair<int, int> > &getYderivDirAndOrder()
        {
            return yDerivDirAndOrder;
        }
        std::vector<std::pair<int, int> > const &getZderivDirAndOrder() const
        {
            return zDerivDirAndOrder;
        }
        std::vector<std::pair<int, int> > &getZderivDirAndOrder()
        {
            return zDerivDirAndOrder;
        }
        Array<T, 3> const &getLocalForce() const
        {
            return localForce;
        }
        Array<T, 3> &getLocalForce()
        {
            return localForce;
        }
        virtual InterpolatedFdOffLatticeInfo3D *clone() const
        {
            return new InterpolatedFdOffLatticeInfo3D(*this);
        }

    private:
        std::vector<Dot3D> wetNodes;
        std::vector<std::vector<std::pair<int, int> > >
            wetNodeSolidUsableDirections;  // stores the directions where there is a wall and valid
                                           // neighbors in its opposite direction
        std::vector<std::vector<int> > wetNodeFluidDirections;  //
        std::vector<std::vector<plint> > wetNodeIds;
        std::vector<std::vector<int> > solidNeighbors;
        std::vector<std::pair<int, int> > xDerivDirAndOrder, yDerivDirAndOrder, zDerivDirAndOrder;
        Array<T, 3> localForce;
    };
};

}  // namespace plb

#endif  // GUO_OFF_LATTICE_MODEL_3D_H
