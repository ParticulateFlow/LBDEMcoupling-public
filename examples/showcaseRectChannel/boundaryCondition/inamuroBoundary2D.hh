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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef INAMURO_BOUNDARY_2D_HH
#define INAMURO_BOUNDARY_2D_HH

#include "boundaryCondition/boundaryInstantiator2D.h"
#include "boundaryCondition/inamuroAnalyticalDynamics.h"
#include "boundaryCondition/inamuroAnalyticalDynamics.hh"
#include "boundaryCondition/inamuroBoundary2D.h"
#include "boundaryCondition/regularizedBoundaryDynamics2D.h"
#include "boundaryCondition/wrappedLocalBoundaryProcessor2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class InamuroBoundaryManager2D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getVelocityBoundaryFunctional();

    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getPressureBoundaryFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getExternalVelocityCornerFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getInternalVelocityCornerFunctional();
};

template <typename T, template <typename U> class Descriptor>
class WrappedInamuroBoundaryManager2D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getVelocityBoundaryFunctional();

    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getPressureBoundaryFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getExternalVelocityCornerFunctional();

    template <int xNormal, int yNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal>
    static BoxProcessingFunctional2D_L<T, Descriptor> *getInternalVelocityCornerFunctional();
};

////////// InamuroBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *InamuroBoundaryManager2D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *InamuroBoundaryManager2D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *InamuroBoundaryManager2D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *InamuroBoundaryManager2D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *InamuroBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreVelocityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *InamuroBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *InamuroBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, xNormal, yNormal>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *InamuroBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return 0;
}

////////// WrappedInamuroBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedInamuroBoundaryManager2D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    bool automaticPrepareCollision = false;
    return new InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedInamuroBoundaryManager2D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedInamuroBoundaryManager2D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    bool automaticPrepareCollision = false;
    return new InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedInamuroBoundaryManager2D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedInamuroBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreVelocityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedInamuroBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new OuterVelocityCornerFunctional2D<T, Descriptor, xNormal, yNormal>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedInamuroBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    bool automaticPrepareCollision = false;
    return new RegularizedVelocityInnerCornerDynamics2D<T, Descriptor, xNormal, yNormal>(
        baseDynamics, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedInamuroBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

////////// Factory function //////////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createInamuroBoundaryCondition2D()
{
    return new BoundaryConditionInstantiator2D<
        T, Descriptor, WrappedInamuroBoundaryManager2D<T, Descriptor> >();
}

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createDynamicsBasedInamuroBoundaryCondition2D()
{
    return new BoundaryConditionInstantiator2D<
        T, Descriptor, InamuroBoundaryManager2D<T, Descriptor> >();
}

}  // namespace plb

#endif  // INAMURO_BOUNDARY_2D_HH
