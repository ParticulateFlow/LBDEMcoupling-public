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
#ifndef GENERALIZED_BOUNDARY_CONDITION_2D_HH
#define GENERALIZED_BOUNDARY_CONDITION_2D_HH

#include "atomicBlock/blockLattice2D.h"
#include "boundaryCondition/boundaryCondition2D.h"
#include "boundaryCondition/boundaryCondition2D.hh"
#include "boundaryCondition/boundaryInstantiator2D.h"
#include "boundaryCondition/generalizedBoundaryCondition2D.h"
#include "core/blockSurface2D.h"
#include "core/plbDebug.h"
#include "generalizedBoundaryDynamics.h"
#include "multiBlock/multiBlockLattice2D.h"

namespace plb {

////////// GeneralizedBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
class GeneralizedBoundaryManager2D {
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
class WrappedGeneralizedBoundaryManager2D {
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

////////// MassConservingGeneralizedBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
class MassConservingGeneralizedBoundaryManager2D {
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
class WrappedMassConservingGeneralizedBoundaryManager2D {
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

////////// GeneralizedBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *GeneralizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoing<Descriptor<T>, direction, orientation>();
    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(baseDynamics, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *GeneralizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *GeneralizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoing<Descriptor<T>, direction, orientation>();
    return new GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *GeneralizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *GeneralizedBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingCorner2D<Descriptor<T>, xNormal, yNormal>();

    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(baseDynamics, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *GeneralizedBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *GeneralizedBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingInternalCorner2D<Descriptor<T>, xNormal, yNormal>();

    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(baseDynamics, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *GeneralizedBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return 0;
}

////////// WrappedGeneralizedBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedGeneralizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryDynamics(
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
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedGeneralizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedGeneralizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoing<Descriptor<T>, direction, orientation>();
    bool automaticPrepareCollision = false;
    return new GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedGeneralizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedGeneralizedBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingCorner2D<Descriptor<T>, xNormal, yNormal>();
    bool automaticPrepareCollision = false;
    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedGeneralizedBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedGeneralizedBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingInternalCorner2D<Descriptor<T>, xNormal, yNormal>();
    bool automaticPrepareCollision = false;

    return new GeneralizedVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor>
    *WrappedGeneralizedBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

////////// MassConservingGeneralizedBoundaryManager2D /////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *MassConservingGeneralizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoing<Descriptor<T>, direction, orientation>();

    std::vector<plint> const &knownIndices =
        indexTemplates::opposite<Descriptor<T> >(missingIndices);
    //         indexTemplates::remainingIndexes<Descriptor<T> >(missingIndices);

    return new GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, knownIndices, missingIndices, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *MassConservingGeneralizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *MassConservingGeneralizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoing<Descriptor<T>, direction, orientation>();

    return new GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor>
    *MassConservingGeneralizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *MassConservingGeneralizedBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingCorner2D<Descriptor<T>, xNormal, yNormal>();
    std::vector<plint> const &knownIndices = indexTemplates::opposite<Descriptor<T> >(
        indexTemplates::subIndexIngoingCorner2D<Descriptor<T>, xNormal, yNormal>());

    std::vector<plint> inGoingIndices;
    inGoingIndices.push_back(
        indexTemplates::findVelocity<Descriptor<T> >(Array<plint, 2>(-xNormal, -yNormal)));

    return new GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, knownIndices, inGoingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor> *
    MassConservingGeneralizedBoundaryManager2D<T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return 0;
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor>
    *MassConservingGeneralizedBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingInternalCorner2D<Descriptor<T>, xNormal, yNormal>();

    return new GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, missingIndices);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor> *
    MassConservingGeneralizedBoundaryManager2D<T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return 0;
}

////////// WrappedMassConservingGeneralizedBoundaryManager2D
////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedMassConservingGeneralizedBoundaryManager2D<T, Descriptor>::getVelocityBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoing<Descriptor<T>, direction, orientation>();

    std::vector<plint> const &knownIndices =
        indexTemplates::opposite<Descriptor<T> >(missingIndices);
    //         indexTemplates::remainingIndexes<Descriptor<T> >(missingIndices);
    bool automaticPrepareCollision = false;

    return new GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, knownIndices, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor> *WrappedMassConservingGeneralizedBoundaryManager2D<
    T, Descriptor>::getVelocityBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoundaryCompositeDynamics<T, Descriptor>
    *WrappedMassConservingGeneralizedBoundaryManager2D<T, Descriptor>::getPressureBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoing<Descriptor<T>, direction, orientation>();
    bool automaticPrepareCollision = false;
    return new GeneralizedDensityBoundaryDynamics<T, Descriptor, direction, orientation>(
        baseDynamics, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int direction, int orientation>
BoxProcessingFunctional2D_L<T, Descriptor> *WrappedMassConservingGeneralizedBoundaryManager2D<
    T, Descriptor>::getPressureBoundaryFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor> *WrappedMassConservingGeneralizedBoundaryManager2D<
    T, Descriptor>::getExternalVelocityCornerDynamics(Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingCorner2D<Descriptor<T>, xNormal, yNormal>();

    std::vector<plint> const &knownIndices = indexTemplates::opposite<Descriptor<T> >(
        indexTemplates::subIndexIngoingCorner2D<Descriptor<T>, xNormal, yNormal>());

    std::vector<plint> inGoingIndices;
    inGoingIndices.push_back(
        indexTemplates::findVelocity<Descriptor<T> >(Array<plint, 2>(-xNormal, -yNormal)));

    bool automaticPrepareCollision = false;
    return new GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, knownIndices, inGoingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor> *WrappedMassConservingGeneralizedBoundaryManager2D<
    T, Descriptor>::getExternalVelocityCornerFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoundaryCompositeDynamics<T, Descriptor> *WrappedMassConservingGeneralizedBoundaryManager2D<
    T, Descriptor>::getInternalVelocityCornerDynamics(Dynamics<T, Descriptor> *baseDynamics)
{
    std::vector<plint> const &missingIndices =
        indexTemplates::subIndexOutgoingInternalCorner2D<Descriptor<T>, xNormal, yNormal>();

    std::vector<plint> const &knownIndices =
        indexTemplates::opposite<Descriptor<T> >(missingIndices);

    bool automaticPrepareCollision = false;

    return new GeneralizedMassConservingVelocityBoundaryDynamics<T, Descriptor>(
        baseDynamics, missingIndices, knownIndices, missingIndices, automaticPrepareCollision);
}

template <typename T, template <typename U> class Descriptor>
template <int xNormal, int yNormal>
BoxProcessingFunctional2D_L<T, Descriptor> *WrappedMassConservingGeneralizedBoundaryManager2D<
    T, Descriptor>::getInternalVelocityCornerFunctional()
{
    return new WrappedLocalBoundaryFunctional2D<T, Descriptor>();
}

////////// Factory functions //////////////////////////////////////////////////

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createGeneralizedBoundaryCondition2D()
{
    return new BoundaryConditionInstantiator2D<
        T, Descriptor, WrappedGeneralizedBoundaryManager2D<T, Descriptor> >;
}

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createDynamicsBasedGeneralizedBoundaryCondition2D()
{
    return new BoundaryConditionInstantiator2D<
        T, Descriptor, GeneralizedBoundaryManager2D<T, Descriptor> >;
}

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createGeneralizedMassConservingBoundaryCondition2D()
{
    return new BoundaryConditionInstantiator2D<
        T, Descriptor, WrappedMassConservingGeneralizedBoundaryManager2D<T, Descriptor> >;
}

template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor>
    *createDynamicsBasedMassConservingGeneralizedBoundaryCondition2D()
{
    return new BoundaryConditionInstantiator2D<
        T, Descriptor, MassConservingGeneralizedBoundaryManager2D<T, Descriptor> >;
}

}  // namespace plb

#endif  // GENERALIZED_BOUNDARY_CONDITION_2D_HH
