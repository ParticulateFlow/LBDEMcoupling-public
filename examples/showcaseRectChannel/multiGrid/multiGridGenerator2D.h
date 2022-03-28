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
 * Various factories that use a multiGridManagement2D -- Header file
 */

#ifndef MULTI_GRID_GENERATOR_2D_H
#define MULTI_GRID_GENERATOR_2D_H

#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiGrid/multiGridLattice2D.h"
#include "multiGrid/multiGridManagement2D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class MultiGridLattice2D;

/// Use the MultiGridManagement2D to generate a vector of lattices that represent the multi grid.
template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlockLattice2D<T, Descriptor> *> generateLattices(
    MultiGridManagement2D management, std::vector<Dynamics<T, Descriptor> *> backgroundDynamics,
    std::vector<BlockCommunicator2D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics, plint envelopeWidth);

template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlockLattice2D<T, Descriptor> *> generateLattices(
    MultiGridManagement2D management, std::vector<Dynamics<T, Descriptor> *> backgroundDynamics,
    plint envelopeWidth);

template <typename T, template <typename U> class Descriptor>
void createInterfaces(
    std::vector<MultiBlockLattice2D<T, Descriptor> *> &multiBlocks,
    MultiGridManagement2D management, bool cubic = false);

/// Instantiate the needed data processors for the coarse grid interface
template <typename T, template <typename U> class Descriptor>
void createCoarseGridInterface(
    plint coarseLevel, Box2D coarseGridInterface,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> &multiBlocks, plint orientation);

/// Instantiate the necessary data processors for the fine grid interface
template <typename T, template <typename U> class Descriptor>
void createFineGridInterface(
    plint coarseLevel, Box2D fineGridInterface,
    std::vector<MultiBlockLattice2D<T, Descriptor> *> &multiBlocks, plint orientation);

/// Use the MultiGridManagement2D to generate a vector of scalar fields with the same distribution
/// as
///  the multi grid
template <typename T>
std::vector<MultiScalarField2D<T> *> generateScalarFields(
    MultiGridManagement2D const &management, std::vector<BlockCommunicator2D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics);

/// Use the MultiGridManagement2D to generate a vector of tensor fields with the same distribution
/// as
///  the multi grid
template <typename T, int nDim>
std::vector<MultiTensorField2D<T, nDim> *> generateTensorFields(
    MultiGridManagement2D const &management, std::vector<BlockCommunicator2D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics);

/// Interface for the class that creates the fine grid interface.
template <typename T, template <typename U> class Descriptor>
class FineGridInterfaceInstantiator {
public:
    virtual ~FineGridInterfaceInstantiator() { }
    virtual FineGridInterfaceInstantiator *clone() = 0;

    virtual void instantiateDataProcessors(
        Box2D whereToInstanciate, MultiBlockLattice2D<T, Descriptor> &coarseLattice,
        MultiBlockLattice2D<T, Descriptor> &fineLattice, plint direction, plint orientation) = 0;

    virtual plint getEnvelopeWidth() = 0;
};

/// Interface for the class that creates the coarse grid interface.
template <typename T, template <typename U> class Descriptor>
class CoarseGridInterfaceInstantiator {
public:
    virtual ~CoarseGridInterfaceInstantiator() { }
    virtual CoarseGridInterfaceInstantiator *clone() = 0;

    virtual void instantiateDataProcessors(
        Box2D whereToInstanciate, MultiBlockLattice2D<T, Descriptor> &coarseLattice,
        MultiBlockLattice2D<T, Descriptor> &fineLattice, plint direction, plint orientation) = 0;
};

/// Create a fine grid interface with linear interpolation.
template <typename T, template <typename U> class Descriptor>
class LinearInterpolationFineGridInterfaceInstantiator :
    public FineGridInterfaceInstantiator<T, Descriptor> {
public:
    virtual LinearInterpolationFineGridInterfaceInstantiator *clone()
    {
        return new LinearInterpolationFineGridInterfaceInstantiator<T, Descriptor>(*this);
    }

    virtual void instantiateDataProcessors(
        Box2D whereToInstanciate, MultiBlockLattice2D<T, Descriptor> &coarseLattice,
        MultiBlockLattice2D<T, Descriptor> &fineLattice, plint direction, plint orientation);

    virtual plint getEnvelopeWidth()
    {
        return 1;
    }
};

/// Create a fine grid interface with cubic interpolation.
template <typename T, template <typename U> class Descriptor>
class CubicInterpolationFineGridInterfaceInstantiator :
    public FineGridInterfaceInstantiator<T, Descriptor> {
public:
    virtual CubicInterpolationFineGridInterfaceInstantiator *clone()
    {
        return new CubicInterpolationFineGridInterfaceInstantiator<T, Descriptor>(*this);
    }

    virtual void instantiateDataProcessors(
        Box2D whereToInstanciate, MultiBlockLattice2D<T, Descriptor> &coarseLattice,
        MultiBlockLattice2D<T, Descriptor> &fineLattice, plint direction, plint orientation);

    virtual plint getEnvelopeWidth()
    {
        return 2;
    }
};

/// Create a coarse grid interface with filtering when copying from fine grid
template <typename T, template <typename U> class Descriptor>
class FilteredCoarseGridInterfaceInstantiator :
    public CoarseGridInterfaceInstantiator<T, Descriptor> {
public:
    virtual FilteredCoarseGridInterfaceInstantiator *clone()
    {
        return new FilteredCoarseGridInterfaceInstantiator<T, Descriptor>(*this);
    }

    virtual void instantiateDataProcessors(
        Box2D whereToInstanciate, MultiBlockLattice2D<T, Descriptor> &coarseLattice,
        MultiBlockLattice2D<T, Descriptor> &fineLattice, plint direction, plint orientation);
};

/// Create a coarse grid interface that performs only a simple copy from the fine grid
template <typename T, template <typename U> class Descriptor>
class UnfilteredCoarseGridInterfaceInstantiator :
    public CoarseGridInterfaceInstantiator<T, Descriptor> {
public:
    virtual UnfilteredCoarseGridInterfaceInstantiator *clone()
    {
        return new UnfilteredCoarseGridInterfaceInstantiator<T, Descriptor>(*this);
    }

    virtual void instantiateDataProcessors(
        Box2D whereToInstanciate, MultiBlockLattice2D<T, Descriptor> &coarseLattice,
        MultiBlockLattice2D<T, Descriptor> &fineLattice, plint direction, plint orientation);
};

/// Factory class to create several different refined grids
template <typename T, template <typename U> class Descriptor>
class MultiGridGenerator2D {
public:
    static std::unique_ptr<MultiGridLattice2D<T, Descriptor> >
        createRefinedLatticeCubicInterpolationNoFiltering(
            MultiGridManagement2D management, Dynamics<T, Descriptor> *backgroundDynamics,
            plint behaviorLevel);

    static std::unique_ptr<MultiGridLattice2D<T, Descriptor> >
        createRefinedLatticeLinearInterpolationNoFiltering(
            MultiGridManagement2D management, Dynamics<T, Descriptor> *backgroundDynamics,
            plint behaviorLevel);

    static std::unique_ptr<MultiGridLattice2D<T, Descriptor> >
        createRefinedLatticeCubicInterpolationFiltering(
            MultiGridManagement2D management, Dynamics<T, Descriptor> *backgroundDynamics,
            plint behaviorLevel);

    static std::unique_ptr<MultiGridLattice2D<T, Descriptor> >
        createRefinedLatticeLinearInterpolationFiltering(
            MultiGridManagement2D management, Dynamics<T, Descriptor> *backgroundDynamics,
            plint behaviorLevel);

private:
    MultiGridGenerator2D();
};

}  // namespace plb

#endif  // MULTI_GRID_GENERATOR_2D_H
