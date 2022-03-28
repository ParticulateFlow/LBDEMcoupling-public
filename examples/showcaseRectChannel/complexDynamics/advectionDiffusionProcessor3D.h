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

#ifndef ADVECTION_DIFFUSION_PROCESSOR_3D_H
#define ADVECTION_DIFFUSION_PROCESSOR_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

/**
 * This class computes the temperature gradient
 * on a convex edge wall in 3D but with a limited number of terms added to the
 * equilibrium distributions (i.e. only the c_i : jNeq term)
 */
template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
class CompleteAdvectionDiffusionEdgeBoundaryFunctional3D :
    public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    enum { direction1 = (plane + 1) % 3, direction2 = (plane + 2) % 3 };

public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual CompleteAdvectionDiffusionEdgeBoundaryFunctional3D<
        T, Descriptor, plane, normal1, normal2>
        *clone() const
    {
        return new CompleteAdvectionDiffusionEdgeBoundaryFunctional3D<
            T, Descriptor, plane, normal1, normal2>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }

public:
    static void processCell(plint iX, plint iY, plint iZ, BlockLattice3D<T, Descriptor> &lattice);

private:
    template <int deriveDirection, int orientation>
    static void interpolateGradients(
        BlockLattice3D<T, Descriptor> const &lattice, T &phiDeriv, plint iX, plint iY, plint iZ);
};

// This data processor uses asymmetric finite differences to compute a gradient.
template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
class CompleteAdvectionDiffusionCornerBoundaryFunctional3D :
    public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);

    virtual CompleteAdvectionDiffusionCornerBoundaryFunctional3D<
        T, Descriptor, xNormal, yNormal, zNormal>
        *clone() const
    {
        return new CompleteAdvectionDiffusionCornerBoundaryFunctional3D<
            T, Descriptor, xNormal, yNormal, zNormal>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;  // lattice
    }
};

// This data processor uses symmetric finite differences to compute a gradient.
// It cannot be applied on any part of the boundary of the global simulation
// domain, except if this boundary is periodic.
template <typename T, template <typename U> class Descriptor>
class SetEffectiveDiffusivity3D : public BoxProcessingFunctional3D_LS<T, Descriptor, T> {
public:
    SetEffectiveDiffusivity3D(T omega0_, T T0_, T cSmago_) :
        omega0(omega0_), invT0((T)1 / T0_), cSmagoSqr(cSmago_ * cSmago_)
    { }
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &rhoBar);
    virtual SetEffectiveDiffusivity3D<T, Descriptor> *clone() const
    {
        return new SetEffectiveDiffusivity3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::dynamicVariables;  // lattice
        modified[1] = modif::nothing;           // rhoBar
    }

private:
    T omega0;
    T invT0;
    T cSmagoSqr;
};

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_PROCESSOR_3D_H
