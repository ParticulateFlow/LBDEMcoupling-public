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

#ifndef GUO_ADV_DIFF_OFF_LATTICE_MODEL_3D_H
#define GUO_ADV_DIFF_OFF_LATTICE_MODEL_3D_H

#include "core/globalDefs.h"
#include "offLattice/offLatticeModel3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class GuoAdvDiffOffLatticeModel3D : public OffLatticeModel3D<T, Array<T, 2> > {
public:
    GuoAdvDiffOffLatticeModel3D(BoundaryShape3D<T, Array<T, 2> > *shape_, int flowType_);
    GuoAdvDiffOffLatticeModel3D(GuoAdvDiffOffLatticeModel3D<T, Descriptor> const &rhs);
    GuoAdvDiffOffLatticeModel3D<T, Descriptor> &operator=(
        GuoAdvDiffOffLatticeModel3D<T, Descriptor> const &rhs);
    virtual GuoAdvDiffOffLatticeModel3D<T, Descriptor> *clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const
    {
        return true;
    }
    virtual void prepareCell(Dot3D const &cellLocation, AtomicContainerBlock3D &container);
    virtual void boundaryCompletion(
        AtomicBlock3D &lattice, AtomicContainerBlock3D &container,
        std::vector<AtomicBlock3D *> const &args);
    virtual ContainerBlockData *generateOffLatticeInfo() const;
    virtual Array<T, 3> getLocalForce(AtomicContainerBlock3D &container) const
    {
        return Array<T, 3>(T(), T(), T());
    }
    void selectSecondOrder(bool flag)
    {
        secondOrderFlag = flag;
    }
    bool usesSecondOrder() const
    {
        return secondOrderFlag;
    }

private:
    void cellCompletion(
        BlockLattice3D<T, Descriptor> &lattice, Dot3D const &guoNode,
        std::vector<std::pair<int, int> > const &dryNodeFluidDirections,
        std::vector<plint> const &dryNodeIds, Dot3D const &absoluteOffset);
    void computeRhoBarJNeq(
        BlockLattice3D<T, Descriptor> const &lattice, Dot3D const &guoNode,
        Dot3D const &fluidDirection, int depth, Array<T, 3> const &wallNode, T delta,
        Array<T, 2> wallData, OffBoundary::Type bdType, Array<T, 3> const &wallNormal, T &rhoBar,
        Array<T, Descriptor<T>::d> &jNeq) const;

private:
    bool secondOrderFlag;

private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class GuoAdvDiffOffLatticeInfo3D : public ContainerBlockData {
    public:
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
        std::vector<std::vector<plint> > const &getDryNodeIds() const
        {
            return dryNodeIds;
        }
        std::vector<std::vector<plint> > &getDryNodeIds()
        {
            return dryNodeIds;
        }
        virtual GuoAdvDiffOffLatticeInfo3D *clone() const
        {
            return new GuoAdvDiffOffLatticeInfo3D(*this);
        }

    private:
        std::vector<Dot3D> dryNodes;
        std::vector<std::vector<std::pair<int, int> > > dryNodeFluidDirections;
        std::vector<std::vector<plint> > dryNodeIds;
    };
};

}  // namespace plb

#endif  // GUO_ADV_DIFF_OFF_LATTICE_MODEL_3D_H
