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

#ifndef FILIPPOVA_HAENEL_OFF_LATTICE_MODEL_3D_H
#define FILIPPOVA_HAENEL_OFF_LATTICE_MODEL_3D_H

#include "core/globalDefs.h"
#include "offLattice/guoOffLatticeModel3D.h"
#include "offLattice/offLatticeModel3D.h"

namespace plb {

/**
 * This class implements the Filippova-Haenel (FH,1998) boundary condition on a BoundaryShape.
 * The BoundaryShape determines whether the points of the discrete lattice are "inside"
 * or "outside" some geometry.
 *
 * It can handle moving boundaries using the momentum correction of ladd (LADD, 1994).
 * The wall velocity is recovered from SurfaceData stored in BoundaryShape3D<T,SurfaceData>*
 *
 * IMPORTANT NOTE: in Palabos versions before June 2020 the FilippovaHaenelLocalModel3D boundary
 * condition refers to the MeiLuoShyy (MLS,1999) variant than now has a independent implementation
 * in offLattice/meiLuoShyyOffLatticeModel3D.h. The name of the class has been changed to have a
 * BREAKING CHANGE in order to encourage the user to choose and check the implementation they need
 * in their application.
 *
 * (FH,1998) O. Filippova and D. Hänel, “Grid Refinement for Lattice-BGK Models,”
 *     Journal of Computational Physics, vol. 147, no. 1, pp. 219–228, Nov. 1998,
 * doi: 10.1006/jcph.1998.6089.
 *
 * (MLS,1999) R. Mei, L.-S. Luo, and W. Shyy, “An Accurate Curved Boundary Treatment in the Lattice
 * Boltzmann Method,” Journal of Computational Physics, vol. 155, no. 2, pp. 307–330, Nov. 1999,
 * doi: 10.1006/jcph.1999.6334.
 *
 * (LADD, 1994) A. J. C. Ladd, “Numerical simulations of particulate suspensions via a discretized
 * Boltzmann equation. Part 1. Theoretical foundation,” Journal of Fluid Mechanics, vol. 271, pp.
 * 285–309, Jul. 1994, doi: 10.1017/S0022112094001771.
 *
 * @tparam T
 * @tparam Descriptor
 */
template <typename T, template <typename U> class Descriptor>
class FilippovaHaenelLocalModel3D : public OffLatticeModel3D<T, Array<T, 3> > {
public:
    FilippovaHaenelLocalModel3D(BoundaryShape3D<T, Array<T, 3> > *shape_, int flowType_);
    virtual FilippovaHaenelLocalModel3D<T, Descriptor> *clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const;
    virtual void prepareCell(Dot3D const &cellLocation, AtomicContainerBlock3D &container);
    virtual void boundaryCompletion(
        AtomicBlock3D &lattice, AtomicContainerBlock3D &container,
        std::vector<AtomicBlock3D *> const &args);

    virtual ContainerBlockData *generateOffLatticeInfo() const;
    virtual Array<T, 3> getLocalForce(AtomicContainerBlock3D &container) const;
    void selectComputeStat(bool flag)
    {
        computeStat = flag;
    }
    bool computesStat() const
    {
        return computeStat;
    }

private:
    void cellCompletion(
        BlockLattice3D<T, Descriptor> &lattice, Dot3D const &guoNode,
        std::vector<int> const &dryNodeFluidDirections, std::vector<plint> const &dryNodeIds,
        Dot3D const &absoluteOffset, Array<T, 3> &localForce,
        std::vector<AtomicBlock3D *> const &args);

private:
    bool computeStat;

private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class OffLatticeInfo3D : public ContainerBlockData {
    public:
        OffLatticeInfo3D() : localForce(Array<T, 3>::zero()) { }
        std::vector<Dot3D> const &getDryNodes() const
        {
            return dryNodes;
        }
        std::vector<Dot3D> &getDryNodes()
        {
            return dryNodes;
        }
        std::vector<std::vector<int> > const &getDryNodeFluidDirections() const
        {
            return dryNodeFluidDirections;
        }
        std::vector<std::vector<int> > &getDryNodeFluidDirections()
        {
            return dryNodeFluidDirections;
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
        virtual OffLatticeInfo3D *clone() const
        {
            return new OffLatticeInfo3D(*this);
        }

    private:
        std::vector<Dot3D> dryNodes;
        std::vector<std::vector<int> > dryNodeFluidDirections;
        std::vector<std::vector<plint> > dryNodeIds;
        Array<T, 3> localForce;
    };
};

}  // namespace plb

#endif  // FILIPPOVA_HAENEL_OFF_LATTICE_MODEL_3D_H
