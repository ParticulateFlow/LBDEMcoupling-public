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

#ifndef BOUZIDI_OFF_LATTICE_MODEL_3D_H
#define BOUZIDI_OFF_LATTICE_MODEL_3D_H

#include "core/globalDefs.h"
#include "offLattice/offLatticeModel3D.h"

namespace plb {

/**
 * This class implements the Bouzidi (BFL,2001) boundary condition on a BoundaryShape.
 * The BoundaryShape determines whether the points of the discrete lattice are voxelFlag::inside
 * or voxelFlag::outside some geometry.
 *
 * It can handle moving boundaries using the momentum correction of ladd (LADD, 1994).
 * The wall velocity is recovered from SurfaceData stored in BoundaryShape3D<T,SurfaceData>*
 *
 * The core of the algorithm is executed by the function cellCompletion().
 *
 * (BFL, 2001) M. Bouzidi, M. Firdaouss, and P. Lallemand, “Momentum transfer of a Boltzmann-lattice
 * fluid with boundaries,” Physics of Fluids, vol. 13, no. 11, pp. 3452–3459, Oct. 2001,
 * doi: 10.1063/1.1399290.
 *
 * (LADD, 1994) A. J. C. Ladd, “Numerical simulations of particulate suspensions via a discretized
 * Boltzmann equation. Part 1. Theoretical foundation,” Journal of Fluid Mechanics, vol. 271, pp.
 * 285–309, Jul. 1994, doi: 10.1017/S0022112094001771.
 *
 * @tparam T
 * @tparam Descriptor
 */
template <typename T, template <typename U> class Descriptor>
class BouzidiOffLatticeModel3D : public OffLatticeModel3D<T, Array<T, 3> > {
public:
    BouzidiOffLatticeModel3D(BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_);
    virtual BouzidiOffLatticeModel3D<T, Descriptor> *clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const;
    virtual void prepareCell(Dot3D const &cellLocation, AtomicContainerBlock3D &container);
    virtual void boundaryCompletion(
        AtomicBlock3D &lattice, AtomicContainerBlock3D &container,
        std::vector<AtomicBlock3D *> const &args);
    void cellCompletion(
        BlockLattice3D<T, Descriptor> &lattice, Dot3D const &boundaryNode,
        std::vector<int> const &solidDirections, std::vector<plint> const &boundaryIds,
        std::vector<bool> const &hasFluidNeighbor, Dot3D const &absoluteOffset,
        Array<T, 3> &localForce, std::vector<AtomicBlock3D *> const &args);
    virtual ContainerBlockData *generateOffLatticeInfo() const;
    virtual Array<T, 3> getLocalForce(AtomicContainerBlock3D &container) const;

private:
    std::vector<T> invAB;

private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class BouzidiOffLatticeInfo3D : public ContainerBlockData {
    public:
        BouzidiOffLatticeInfo3D() : localForce(Array<T, 3>::zero()) { }
        std::vector<Dot3D> const &getBoundaryNodes() const
        {
            return boundaryNodes;
        }
        std::vector<Dot3D> &getBoundaryNodes()
        {
            return boundaryNodes;
        }
        std::vector<std::vector<int> > const &getSolidDirections() const
        {
            return solidDirections;
        }
        std::vector<std::vector<int> > &getSolidDirections()
        {
            return solidDirections;
        }
        std::vector<std::vector<plint> > const &getBoundaryIds() const
        {
            return boundaryIds;
        }
        std::vector<std::vector<plint> > &getBoundaryIds()
        {
            return boundaryIds;
        }
        std::vector<std::vector<bool> > const &getHasFluidNeighbor() const
        {
            return hasFluidNeighbor;
        }
        std::vector<std::vector<bool> > &getHasFluidNeighbor()
        {
            return hasFluidNeighbor;
        }
        Array<T, 3> const &getLocalForce() const
        {
            return localForce;
        }
        Array<T, 3> &getLocalForce()
        {
            return localForce;
        }
        virtual BouzidiOffLatticeInfo3D *clone() const
        {
            return new BouzidiOffLatticeInfo3D(*this);
        }

    private:
        std::vector<Dot3D> boundaryNodes;
        std::vector<std::vector<int> > solidDirections;
        std::vector<std::vector<plint> > boundaryIds;
        std::vector<std::vector<bool> > hasFluidNeighbor;
        Array<T, 3> localForce;
    };
};

}  // namespace plb

#endif  // BOUZIDI_OFF_LATTICE_MODEL_3D_H
