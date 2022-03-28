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
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_HH
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_HH

#include "complexDynamics/advectionDiffusionBoundaries.h"
#include "complexDynamics/advectionDiffusionBoundaries.hh"
#include "complexDynamics/advectionDiffusionBoundaryCondition2D.h"
#include "complexDynamics/advectionDiffusionBoundaryInstantiator2D.h"
#include "complexDynamics/advectionDiffusionProcessor2D.h"
#include "complexDynamics/advectionDiffusionProcessor2D.hh"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void OnLatticeAdvectionDiffusionBoundaryCondition2D<
    T, Descriptor>::setTemperatureConditionOnBlockBoundaries(BlockLattice2D<T, Descriptor> &lattice)
{
    plint nx = lattice.getNx();
    plint ny = lattice.getNy();

    addTemperatureBoundary0N(Box2D(0, 0, 1, ny - 2), lattice);
    addTemperatureBoundary0P(Box2D(nx - 1, nx - 1, 1, ny - 2), lattice);
    addTemperatureBoundary1N(Box2D(1, nx - 2, 0, 0), lattice);
    addTemperatureBoundary1P(Box2D(1, nx - 2, ny - 1, ny - 1), lattice);

    addTemperatureCornerNN(0, 0, lattice);
    addTemperatureCornerNP(0, ny - 1, lattice);
    addTemperatureCornerPN(nx - 1, 0, lattice);
    addTemperatureCornerPP(nx - 1, ny - 1, lattice);
}

template <typename T, template <typename U> class Descriptor>
void OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor>::
    setTemperatureConditionOnBlockBoundaries(MultiBlockLattice2D<T, Descriptor> &lattice)
{
    plint nx = lattice.getNx();
    plint ny = lattice.getNy();

    addTemperatureBoundary0N(Box2D(0, 0, 1, ny - 2), lattice);
    addTemperatureBoundary0P(Box2D(nx - 1, nx - 1, 1, ny - 2), lattice);
    addTemperatureBoundary1N(Box2D(1, nx - 2, 0, 0), lattice);
    addTemperatureBoundary1P(Box2D(1, nx - 2, ny - 1, ny - 1), lattice);

    addTemperatureCornerNN(0, 0, lattice);
    addTemperatureCornerNP(0, ny - 1, lattice);
    addTemperatureCornerPN(nx - 1, 0, lattice);
    addTemperatureCornerPP(nx - 1, ny - 1, lattice);
}

//================ Zou-He like advectionDiffusionBoundaryManager2D ==========//

template <typename T, template <typename U> class Descriptor>
class AdvectionDiffusionBoundaryManager2D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getAdvectionDiffusionBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getAdvectionDiffusionBoundaryProcessor();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getAdvectionDiffusionCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getAdvectionDiffusionCornerProcessor();
};

//================ Regularized advectionDiffusionBoundaryManager2D ==========//

template <typename T, template <typename U> class Descriptor>
class RegularizedAdvectionDiffusionBoundaryManager2D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getAdvectionDiffusionBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getAdvectionDiffusionBoundaryProcessor();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getAdvectionDiffusionCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getAdvectionDiffusionCornerProcessor();
};

//================ Complete Regularized advectionDiffusionBoundaryManager2D ==========//

template <typename T, template <typename U> class Descriptor>
class RegularizedCompleteAdvectionDiffusionBoundaryManager2D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getAdvectionDiffusionBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getAdvectionDiffusionBoundaryProcessor();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getAdvectionDiffusionCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getAdvectionDiffusionCornerProcessor();
};

////////// Zou-He AdvectionDiffusionBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *AdvectionDiffusionBoundaryManager2D<T, Descriptor>::getAdvectionDiffusionBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *AdvectionDiffusionBoundaryManager2D<T, Descriptor>::getAdvectionDiffusionBoundaryProcessor()
{
    return 0;
}

//==================  Corners ================================

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *AdvectionDiffusionBoundaryManager2D<T, Descriptor>::getAdvectionDiffusionCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *AdvectionDiffusionBoundaryManager2D<T, Descriptor>::getAdvectionDiffusionCornerProcessor()
{
    return 0;
}

////////// Regularized AdvectionDiffusionBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor> *RegularizedAdvectionDiffusionBoundaryManager2D<
    T, Descriptor>::getAdvectionDiffusionBoundaryDynamics(Dynamics<T, Descriptor> *baseDynamics)
{
    return new RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor> *RegularizedAdvectionDiffusionBoundaryManager2D<
    T, Descriptor>::getAdvectionDiffusionBoundaryProcessor()
{
    return 0;
}

//==================  Corners ================================

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor> *RegularizedAdvectionDiffusionBoundaryManager2D<
    T, Descriptor>::getAdvectionDiffusionCornerDynamics(Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreDensityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor> *RegularizedAdvectionDiffusionBoundaryManager2D<
    T, Descriptor>::getAdvectionDiffusionCornerProcessor()
{
    return new AdvectionDiffusionCornerBoundaryFunctional2D<T, Descriptor, xNormal, yNormal>();
}

////////// Regularized CompleteAdvectionDiffusionBoundaryManager2D
////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor> *RegularizedCompleteAdvectionDiffusionBoundaryManager2D<
    T, Descriptor>::getAdvectionDiffusionBoundaryDynamics(Dynamics<T, Descriptor> *baseDynamics)
{
    return new CompleteRegularizedAdvectionDiffusionBoundaryDynamics<
        T, Descriptor, direction, orientation>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor> *RegularizedCompleteAdvectionDiffusionBoundaryManager2D<
    T, Descriptor>::getAdvectionDiffusionBoundaryProcessor()
{
    return 0;
}

//==================  Corners ================================

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor> *RegularizedCompleteAdvectionDiffusionBoundaryManager2D<
    T, Descriptor>::getAdvectionDiffusionCornerDynamics(Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreDensityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor> *RegularizedCompleteAdvectionDiffusionBoundaryManager2D<
    T, Descriptor>::getAdvectionDiffusionCornerProcessor()
{
    return new CompleteAdvectionDiffusionCornerBoundaryFunctional2D<
        T, Descriptor, xNormal, yNormal>();
}

////////// Factory functions //////////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor>
    *createLocalAdvectionDiffusionBoundaryCondition2D()
{
    return new AdvectionDiffusionBoundaryConditionInstantiator2D<
        T, Descriptor, AdvectionDiffusionBoundaryManager2D<T, Descriptor> >();
}

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor>
    *createLocalRegularizedAdvectionDiffusionBoundaryCondition2D()
{
    return new AdvectionDiffusionBoundaryConditionInstantiator2D<
        T, Descriptor, RegularizedAdvectionDiffusionBoundaryManager2D<T, Descriptor> >();
}

template <typename T, template <typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Descriptor>
    *createLocalRegularizedCompleteAdvectionDiffusionBoundaryCondition2D()
{
    return new AdvectionDiffusionBoundaryConditionInstantiator2D<
        T, Descriptor, RegularizedCompleteAdvectionDiffusionBoundaryManager2D<T, Descriptor> >();
}

}  // namespace plb

#endif
