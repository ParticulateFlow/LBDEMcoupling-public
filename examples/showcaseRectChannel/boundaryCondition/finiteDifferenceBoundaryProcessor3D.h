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

#ifndef FINITE_DIFFERENCE_BOUNDARY_PROCESSOR_3D_H
#define FINITE_DIFFERENCE_BOUNDARY_PROCESSOR_3D_H

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/globalDefs.h"

namespace plb {

/**
 * This class computes the Skordos BC
 * on a plane wall in 3D but with a limited number of terms added to the
 * equilibrium distributions (i.e. only the Q_i : Pi term)
 */
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class PlaneFdBoundaryFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual PlaneFdBoundaryFunctional3D<T, Descriptor, direction, orientation> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
    virtual int getStaticId() const
    {
        return staticId;
    }
    static const int staticId;

public:
    static void processCell(plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice);

private:
    template <int deriveDirection>
    static void interpolateGradients(
        BlockLattice3D<T, Descriptor> const &block, Array<T, Descriptor<T>::d> &velDeriv, plint iX,
        plint iY, plint iZ);
};

/**
 * This class computes the Skordos BC
 * on a convex edge wall in 3D but with a limited number of terms added to the
 * equilibrium distributions (i.e. only the Q_i : Pi term)
 */
template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
class OuterVelocityEdgeFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    enum { direction1 = (plane + 1) % 3, direction2 = (plane + 2) % 3 };

public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
    virtual int getStaticId() const
    {
        return staticId;
    }
    static const int staticId;

public:
    static void processCell(plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice);

private:
    static T getNeighborRho(
        plint x, plint y, plint z, plint step1, plint step2,
        BlockLattice3D<T, Descriptor> const &lattice);
    template <int deriveDirection, int orientation>
    static void interpolateGradients(
        BlockLattice3D<T, Descriptor> const &lattice, Array<T, Descriptor<T>::d> &velDeriv,
        plint iX, plint iY, plint iZ);
};

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
class OuterVelocityCornerFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual OuterVelocityCornerFunctional3D<T, Descriptor, xNormal, yNormal, zNormal> *clone()
        const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
    virtual int getStaticId() const
    {
        return staticId;
    }

public:
    static void processCell(plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice);
    static const int staticId;
};

}  // namespace plb

#endif  // FINITE_DIFFERENCE_BOUNDARY_PROCESSOR_3D_H
