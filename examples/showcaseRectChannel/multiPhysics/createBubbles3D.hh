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

#ifndef CREATE_BUBBLES_3D_HH
#define CREATE_BUBBLES_3D_HH

#include <limits>

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "multiPhysics/createBubbles3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"
#include "offLattice/makeSparse3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void punchSphere(
    FreeSurfaceFields3D<T, Descriptor> &fields, Array<T, 3> const &center, T radius, T rhoEmpty,
    T rho0, Dynamics<T, Descriptor> &dynamics)
{
    applyProcessingFunctional(
        new PunchSphere3D<T, Descriptor>(center, radius, rho0), fields.lattice.getBoundingBox(),
        fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), fields.lattice.getBoundingBox(),
        fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoEmpty),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(
            dynamics.clone(), Array<T, 3>((T)0., (T)0., (T)0.)),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoEmpty),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>(),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceComputeStatistics3D<T, Descriptor>, fields.lattice.getBoundingBox(),
        fields.freeSurfaceArgs);
}

template <typename T, template <typename U> class Descriptor>
void analyticalPunchSphere(
    FreeSurfaceFields3D<T, Descriptor> &fields, Array<T, 3> const &center, T radius, T rhoEmpty,
    T rho0, plint subDivision, Dynamics<T, Descriptor> &dynamics)
{
    applyProcessingFunctional(
        new AnalyticalPunchSphere3D<T, Descriptor>(center, radius, rho0, subDivision),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), fields.lattice.getBoundingBox(),
        fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoEmpty),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(
            dynamics.clone(), Array<T, 3>((T)0., (T)0., (T)0.)),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoEmpty),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>(),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceComputeStatistics3D<T, Descriptor>, fields.lattice.getBoundingBox(),
        fields.freeSurfaceArgs);
}

template <typename T, template <typename U> class Descriptor>
T computeAverageSphereDensity(
    FreeSurfaceFields3D<T, Descriptor> &fields, Array<T, 3> const &center, T radius)
{
    CalculateAverageSphereDensity3D<T, Descriptor> functional(center, radius);
    applyProcessingFunctional(functional, fields.lattice.getBoundingBox(), fields.freeSurfaceArgs);
    return functional.getAverageDensity();
}

template <typename T, template <typename U> class Descriptor>
void punchSphere(
    FreeSurfaceSetup<T, Descriptor> &setup, Array<T, 3> const &center, T radius, T rhoEmpty, T rho0,
    Dynamics<T, Descriptor> &dynamics)
{
    applyProcessingFunctional(
        new PunchSphere3D<T, Descriptor>(center, radius, rho0), setup.getGroup().getBoundingBox(),
        setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), setup.getGroup().getBoundingBox(),
        setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoEmpty),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(
            dynamics.clone(), Array<T, 3>((T)0., (T)0., (T)0.)),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoEmpty),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>(),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceComputeStatistics3D<T, Descriptor>, setup.getGroup().getBoundingBox(),
        setup.freeSurfaceArgs);
}

template <typename T, template <typename U> class Descriptor>
void analyticalPunchSphere(
    FreeSurfaceSetup<T, Descriptor> &setup, Array<T, 3> const &center, T radius, T rhoEmpty, T rho0,
    plint subDivision, Dynamics<T, Descriptor> &dynamics)
{
    applyProcessingFunctional(
        new AnalyticalPunchSphere3D<T, Descriptor>(center, radius, rho0, subDivision),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), setup.getGroup().getBoundingBox(),
        setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoEmpty),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(
            dynamics.clone(), Array<T, 3>((T)0., (T)0., (T)0.)),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoEmpty),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>(),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceComputeStatistics3D<T, Descriptor>, setup.getGroup().getBoundingBox(),
        setup.freeSurfaceArgs);
}

template <typename T, template <typename U> class Descriptor>
T computeAverageSphereDensity(
    FreeSurfaceSetup<T, Descriptor> &setup, Array<T, 3> const &center, T radius)
{
    CalculateAverageSphereDensity3D<T, Descriptor> functional(center, radius);
    applyProcessingFunctional(functional, setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs);
    return functional.getAverageDensity();
}

template <typename T, template <typename U> class Descriptor>
void punchSphere(
    FreeSurfaceWrapper<T, Descriptor> &wrapper, Array<T, 3> const &center, T radius, T rhoEmpty,
    T rho0, Dynamics<T, Descriptor> &dynamics)
{
    applyProcessingFunctional(
        new PunchSphere3D<T, Descriptor>(center, radius, rho0), wrapper.getGroup().getBoundingBox(),
        wrapper.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoEmpty),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(
            dynamics.clone(), Array<T, 3>((T)0., (T)0., (T)0.)),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoEmpty),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceEqualMassExcessReDistributionAndComputationOfLostMass3D<T, Descriptor>(
            wrapper.getReductionData()),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    {
        std::vector<MultiBlock3D *> args;
        args.push_back(&wrapper.getFlag());
        args.push_back(&wrapper.getMass());

        applyProcessingFunctional(
            new FreeSurfaceComputeReductionsPerProcess3D<T>(wrapper.getReductionData()),
            wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);
    }

    {
        std::vector<MultiBlock3D *> args;
        args.push_back(&wrapper.getLattice());

#ifdef PLB_MPI_PARALLEL
        applyProcessingFunctional(
            new FreeSurfaceComputeReductions3D(
                wrapper.getReductionData(), wrapper.getReductionCommunicator()),
            wrapper.getGroup().getBoundingBox(), args);
#else
        applyProcessingFunctional(
            new FreeSurfaceComputeSerialReductions3D(wrapper.getReductionData()),
            wrapper.getGroup().getBoundingBox(), args);
#endif

        applyProcessingFunctional(
            new FreeSurfaceResetReductionData3D(wrapper.getReductionData()),
            wrapper.getGroup().getBoundingBox(), args);
    }
}

template <typename T, template <typename U> class Descriptor>
void analyticalPunchSphere(
    FreeSurfaceWrapper<T, Descriptor> &wrapper, Array<T, 3> const &center, T radius, T rhoEmpty,
    T rho0, plint subDivision, Dynamics<T, Descriptor> &dynamics)
{
    applyProcessingFunctional(
        new AnalyticalPunchSphere3D<T, Descriptor>(center, radius, rho0, subDivision),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoEmpty),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(
            dynamics.clone(), Array<T, 3>((T)0., (T)0., (T)0.)),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoEmpty),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    applyProcessingFunctional(
        new FreeSurfaceEqualMassExcessReDistributionAndComputationOfLostMass3D<T, Descriptor>(
            wrapper.getReductionData()),
        wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);

    {
        std::vector<MultiBlock3D *> args;
        args.push_back(&wrapper.getFlag());
        args.push_back(&wrapper.getMass());

        applyProcessingFunctional(
            new FreeSurfaceComputeReductionsPerProcess3D<T>(wrapper.getReductionData()),
            wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);
    }

    {
        std::vector<MultiBlock3D *> args;
        args.push_back(&wrapper.getLattice());

#ifdef PLB_MPI_PARALLEL
        applyProcessingFunctional(
            new FreeSurfaceComputeReductions3D(
                wrapper.getReductionData(), wrapper.getReductionCommunicator()),
            wrapper.getGroup().getBoundingBox(), args);
#else
        applyProcessingFunctional(
            new FreeSurfaceComputeSerialReductions3D(wrapper.getReductionData()),
            wrapper.getGroup().getBoundingBox(), args);
#endif

        applyProcessingFunctional(
            new FreeSurfaceResetReductionData3D(wrapper.getReductionData()),
            wrapper.getGroup().getBoundingBox(), args);
    }
}

template <typename T, template <typename U> class Descriptor>
T computeAverageSphereDensity(
    FreeSurfaceWrapper<T, Descriptor> &wrapper, Array<T, 3> const &center, T radius)
{
    CalculateAverageSphereDensity3D<T, Descriptor> functional(center, radius);
    applyProcessingFunctional(
        functional, wrapper.getGroup().getBoundingBox(), wrapper.freeSurfaceArgs);
    return functional.getAverageDensity();
}

}  // namespace plb

#endif  // CREATE_BUBBLES_3D_HH
