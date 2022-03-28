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

#ifndef ZOU_HE_BOUNDARY_3D_HH
#define ZOU_HE_BOUNDARY_3D_HH

#include "boundaryCondition/boundaryInstantiator3D.h"
#include "boundaryCondition/regularizedBoundaryDynamics3D.h"
#include "boundaryCondition/wrappedLocalBoundaryProcessor3D.h"
#include "boundaryCondition/zouHeBoundary3D.h"
#include "boundaryCondition/zouHeDynamics.h"
#include "boundaryCondition/zouHeDynamics.hh"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class ZouHeBoundaryManager3D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getVelocityBoundaryFunctional();

    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getPressureBoundaryFunctional();

    template <int plane, int normal1, int normal2>
    static BoundaryCompositeDynamics<T, Descriptor> *getExternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int plane, int normal1, int normal2>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getExternalVelocityEdgeFunctional();

    template <int plane, int normal1, int normal2>
    static BoundaryCompositeDynamics<T, Descriptor> *getInternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int plane, int normal1, int normal2>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getInternalVelocityEdgeFunctional();

    template <int xNormal, int yNormal, int zNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal, int zNormal>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getExternalVelocityCornerFunctional();

    template <int xNormal, int yNormal, int zNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal, int zNormal>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getInternalVelocityCornerFunctional();
};

template <typename T, template <typename U> class Descriptor>
class WrappedZouHeBoundaryManager3D {
public:
    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getVelocityBoundaryFunctional();

    template <int direction, int orientation>
    static BoundaryCompositeDynamics<T, Descriptor> *getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int direction, int orientation>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getPressureBoundaryFunctional();

    template <int plane, int normal1, int normal2>
    static BoundaryCompositeDynamics<T, Descriptor> *getExternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int plane, int normal1, int normal2>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getExternalVelocityEdgeFunctional();

    template <int plane, int normal1, int normal2>
    static BoundaryCompositeDynamics<T, Descriptor> *getInternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int plane, int normal1, int normal2>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getInternalVelocityEdgeFunctional();

    template <int xNormal, int yNormal, int zNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal, int zNormal>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getExternalVelocityCornerFunctional();

    template <int xNormal, int yNormal, int zNormal>
    static BoundaryCompositeDynamics<T, Descriptor> *getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics);
    template <int xNormal, int yNormal, int zNormal>
    static BoxProcessingFunctional3D_L<T, Descriptor> *getInternalVelocityCornerFunctional();
};

////////// ZouHeBoundaryManager3D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new ZouHeVelocityDynamics<T, Descriptor, direction, orientation>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new ZouHePressureDynamics<T, Descriptor, direction, orientation>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getExternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreVelocityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getExternalVelocityEdgeFunctional()
{
    return new OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>();
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getInternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getInternalVelocityEdgeFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreVelocityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new OuterVelocityCornerFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>(
        baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor>
    *ZouHeBoundaryManager3D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return 0;
}

////////// WrappedZouHeBoundaryManager3D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    bool automaticPrepareCollision = false;
    return new ZouHeVelocityDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional3D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    bool automaticPrepareCollision = false;
    return new ZouHePressureDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional3D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getExternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreVelocityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getExternalVelocityEdgeFunctional()
{
    return new OuterVelocityEdgeFunctional3D<T, Descriptor, plane, normal1, normal2>();
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getInternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    bool automaticPrepareCollision = false;
    return new RegularizedVelocityInnerEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>(
        baseDynamics, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getInternalVelocityEdgeFunctional()
{
    return new WrappedLocalBoundaryFunctional3D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    return new StoreVelocityDynamics<T, Descriptor>(baseDynamics);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new OuterVelocityCornerFunctional3D<T, Descriptor, xNormal, yNormal, zNormal>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    bool automaticPrepareCollision = false;
    return new RegularizedVelocityInnerCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>(
        baseDynamics, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedZouHeBoundaryManager3D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return new WrappedLocalBoundaryFunctional3D<T, Descriptor>();
}

////////// Factory functions //////////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T, Descriptor> *createZouHeBoundaryCondition3D()
{
    return new BoundaryConditionInstantiator3D<
        T, Descriptor, WrappedZouHeBoundaryManager3D<T, Descriptor> >();
}

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T, Descriptor> *createDynamicsBasedZouHeBoundaryCondition3D()
{
    return new BoundaryConditionInstantiator3D<
        T, Descriptor, ZouHeBoundaryManager3D<T, Descriptor> >();
}

}  // namespace plb

#endif  // ZOU_HE_BOUNDARY_3D_HH
