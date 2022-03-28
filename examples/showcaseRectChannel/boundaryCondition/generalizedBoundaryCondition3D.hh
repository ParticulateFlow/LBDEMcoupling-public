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
 * A helper for initialising 3D boundaries -- generic implementation.
 */
#ifndef GENERALIZED_BOUNDARY_CONDITION_3D_HH
#define GENERALIZED_BOUNDARY_CONDITION_3D_HH

#include "boundaryCondition/boundaryCondition3D.h"
#include "boundaryCondition/boundaryInstantiator3D.h"
#include "boundaryCondition/equilibriumBoundaryDynamics.h"
#include "boundaryCondition/regularizedBoundaryDynamics3D.h"
#include "core/blockSurface3D.h"
#include "core/plbDebug.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class GeneralizedBoundaryManager3D {
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
class WrappedGeneralizedBoundaryManager3D {
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

////////// GeneralizedBoundaryManager3D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoing<Descriptor<T>, direction, orientation>();

    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(baseDynamics, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    //     TODO implement pressure BC
    PLB_ASSERT(false);
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getExternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingExternalEdge3D<Descriptor<T>, plane, normal1, normal2>();
    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(baseDynamics, missingIndices);
    ;
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getExternalVelocityEdgeFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getInternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingInternalEdge3D<Descriptor<T>, plane, normal1, normal2>();
    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(baseDynamics, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getInternalVelocityEdgeFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices = indexTemplates::subIndexOutgoingExternalCorner3D<
        Descriptor<T>, xNormal, yNormal, zNormal>();

    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(baseDynamics, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices = indexTemplates::subIndexOutgoingInternalCorner3D<
        Descriptor<T>, xNormal, yNormal, zNormal>();
    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(baseDynamics, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor>
    *GeneralizedBoundaryManager3D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return 0;
}

////////// WrappedGeneralizedBoundaryManager3D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoing<Descriptor<T>, direction, orientation>();
    bool automaticPrepareCollision = false;

    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional3D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    //     TODO implement pressure BC
    PLB_ASSERT(false);
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getExternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingExternalEdge3D<Descriptor<T>, plane, normal1, normal2>();
    bool automaticPrepareCollision = false;
    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getExternalVelocityEdgeFunctional()
{
    return new WrappedLocalBoundaryFunctional3D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getInternalVelocityEdgeDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingInternalEdge3D<Descriptor<T>, plane, normal1, normal2>();
    bool automaticPrepareCollision = false;
    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int plane, int normal1, int normal2>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getInternalVelocityEdgeFunctional()
{
    return new WrappedLocalBoundaryFunctional3D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices = indexTemplates::subIndexOutgoingExternalCorner3D<
        Descriptor<T>, xNormal, yNormal, zNormal>();
    bool automaticPrepareCollision = false;

    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new WrappedLocalBoundaryFunctional3D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices = indexTemplates::subIndexOutgoingInternalCorner3D<
        Descriptor<T>, xNormal, yNormal, zNormal>();
    bool automaticPrepareCollision = false;
    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal, int zNormal>
BoxProcessingFunctional3D_L<T, Descriptor>
    *WrappedGeneralizedBoundaryManager3D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return new WrappedLocalBoundaryFunctional3D<T, Descriptor>();
}

////////// Factory functions //////////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T, Descriptor> *createGeneralizedBoundaryCondition3D()
{
    return new BoundaryConditionInstantiator3D<
        T, Descriptor, WrappedGeneralizedBoundaryManager3D<T, Descriptor> >();
}

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T, Descriptor> *createDynamicsBasedGeneralizedBoundaryCondition3D()
{
    return new BoundaryConditionInstantiator3D<
        T, Descriptor, GeneralizedBoundaryManager3D<T, Descriptor> >();
}

}  // namespace plb

#endif  // GENERALIZED_BOUNDARY_CONDITION_3D_HH
