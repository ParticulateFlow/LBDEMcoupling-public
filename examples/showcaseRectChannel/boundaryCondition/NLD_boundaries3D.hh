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
 * A framework for non-local boundary conditions instantiated in terms of
 * local dynamics objects -- generic implementation.
 */
#ifndef NLD_BOUNDARIES_3D_HH
#define NLD_BOUNDARIES_3D_HH

#include "boundaryCondition/NLD_boundaries3D.h"
#include "core/blockSurface3D.h"
#include "core/globalDefs.h"
#include "core/nonLocalDynamics3D.h"
#include "core/processorIdentifiers3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"

namespace plb {

////////  ExecutePlaneNLD_3D ///////////////////////////////////

template <typename T, template <typename U> class Descriptor>
const int ExecutePlaneNLD_3D<T, Descriptor>::staticId =
    meta::registerProcessor3D<ExecutePlaneNLD_3D<T, Descriptor>, T, Descriptor>(
        std::string("ExecutePlaneNLD_3D"));

template <typename T, template <typename U> class Descriptor>
ExecutePlaneNLD_3D<T, Descriptor>::ExecutePlaneNLD_3D() : direction(0), orientation(0)
{ }

template <typename T, template <typename U> class Descriptor>
ExecutePlaneNLD_3D<T, Descriptor>::ExecutePlaneNLD_3D(int direction_, int orientation_) :
    direction(direction_), orientation(orientation_)
{ }

template <typename T, template <typename U> class Descriptor>
void ExecutePlaneNLD_3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Dynamics<T, Descriptor> *originalDynamics = &lattice.get(iX, iY, iZ).getDynamics();
                if (originalDynamics->isNonLocal()) {
                    NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics =
                        (NonLocalBoundaryDynamics3D<T, Descriptor> *)(originalDynamics);
                    dynamics->planeBoundaryCompletion(direction, orientation, iX, iY, iZ, lattice);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ExecutePlaneNLD_3D<T, Descriptor> *ExecutePlaneNLD_3D<T, Descriptor>::clone() const
{
    return new ExecutePlaneNLD_3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExecutePlaneNLD_3D<T, Descriptor>::serialize(std::string &data) const
{
    std::string newData = util::val2str(direction, orientation);
    data += newData + " ";
}

template <typename T, template <typename U> class Descriptor>
void ExecutePlaneNLD_3D<T, Descriptor>::unserialize(std::string &data)
{
    data = util::consume(data, direction, orientation);
}

////////  ExecuteEdgeNLD_3D ///////////////////////////////////

template <typename T, template <typename U> class Descriptor>
const int ExecuteEdgeNLD_3D<T, Descriptor>::staticId =
    meta::registerProcessor3D<ExecuteEdgeNLD_3D<T, Descriptor>, T, Descriptor>(
        std::string("ExecuteEdgeNLD_3D"));

template <typename T, template <typename U> class Descriptor>
ExecuteEdgeNLD_3D<T, Descriptor>::ExecuteEdgeNLD_3D() : plane(0), normal1(0), normal2(0)
{ }

template <typename T, template <typename U> class Descriptor>
ExecuteEdgeNLD_3D<T, Descriptor>::ExecuteEdgeNLD_3D(int plane_, int normal1_, int normal2_) :
    plane(plane_), normal1(normal1_), normal2(normal2_)
{ }

template <typename T, template <typename U> class Descriptor>
void ExecuteEdgeNLD_3D<T, Descriptor>::process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics =
                    dynamic_cast<NonLocalBoundaryDynamics3D<T, Descriptor> *>(
                        &lattice.get(iX, iY, iZ).getDynamics());
                if (dynamics) {
                    dynamics->edgeBoundaryCompletion(plane, normal1, normal2, iX, iY, iZ, lattice);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ExecuteEdgeNLD_3D<T, Descriptor> *ExecuteEdgeNLD_3D<T, Descriptor>::clone() const
{
    return new ExecuteEdgeNLD_3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExecuteEdgeNLD_3D<T, Descriptor>::serialize(std::string &data) const
{
    std::string newData = util::val2str(plane, normal1, normal2);
    data += newData + " ";
}

template <typename T, template <typename U> class Descriptor>
void ExecuteEdgeNLD_3D<T, Descriptor>::unserialize(std::string &data)
{
    data = util::consume(data, plane, normal1, normal2);
}

////////  ExecuteCornerNLD_3D ///////////////////////////////////

template <typename T, template <typename U> class Descriptor>
const int ExecuteCornerNLD_3D<T, Descriptor>::staticId =
    meta::registerProcessor3D<ExecuteCornerNLD_3D<T, Descriptor>, T, Descriptor>(
        std::string("ExecuteCornerNLD_3D"));

template <typename T, template <typename U> class Descriptor>
ExecuteCornerNLD_3D<T, Descriptor>::ExecuteCornerNLD_3D() : xNormal(0), yNormal(0), zNormal(0)
{ }

template <typename T, template <typename U> class Descriptor>
ExecuteCornerNLD_3D<T, Descriptor>::ExecuteCornerNLD_3D(int xNormal_, int yNormal_, int zNormal_) :
    xNormal(xNormal_), yNormal(yNormal_), zNormal(zNormal_)
{ }

template <typename T, template <typename U> class Descriptor>
void ExecuteCornerNLD_3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics =
                    dynamic_cast<NonLocalBoundaryDynamics3D<T, Descriptor> *>(
                        &lattice.get(iX, iY, iZ).getDynamics());
                if (dynamics) {
                    dynamics->cornerBoundaryCompletion(
                        xNormal, yNormal, zNormal, iX, iY, iZ, lattice);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ExecuteCornerNLD_3D<T, Descriptor> *ExecuteCornerNLD_3D<T, Descriptor>::clone() const
{
    return new ExecuteCornerNLD_3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ExecuteCornerNLD_3D<T, Descriptor>::serialize(std::string &data) const
{
    std::string newData = util::val2str(xNormal, yNormal, zNormal);
    data += newData + " ";
}

template <typename T, template <typename U> class Descriptor>
void ExecuteCornerNLD_3D<T, Descriptor>::unserialize(std::string &data)
{
    data = util::consume(data, xNormal, yNormal, zNormal);
}

////////  ExecuteNonLocalDynamics3D ///////////////////////////////////

template <typename T, template <typename U> class Descriptor>
const int ExecuteNonLocalDynamics3D<T, Descriptor>::staticId =
    meta::registerProcessor3D<ExecuteNonLocalDynamics3D<T, Descriptor>, T, Descriptor>(
        std::string("ExecuteNonLocal3D"));

template <typename T, template <typename U> class Descriptor>
void ExecuteNonLocalDynamics3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Dynamics<T, Descriptor> &dynamics = lattice.get(iX, iY, iZ).getDynamics();
                if (dynamics.isNonLocal()) {
                    NonLocalDynamics3D<T, Descriptor> &nonLocalDynamics(
                        dynamic_cast<NonLocalDynamics3D<T, Descriptor> &>(dynamics));
                    nonLocalDynamics.nonLocalAction(iX, iY, iZ, lattice);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ExecuteNonLocalDynamics3D<T, Descriptor> *ExecuteNonLocalDynamics3D<T, Descriptor>::clone() const
{
    return new ExecuteNonLocalDynamics3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void setFluidNLDboundaryDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, boundary::BcType bcType)
{
    setFluidNLDboundaryDynamics(lattice, lattice.getBoundingBox(), domain, bcType);
}

template <typename T, template <typename U> class Descriptor>
void setFluidNLDboundaryDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D bbox, Box3D domain, boundary::BcType bcType)
{
    Box3D domain1(bbox.x0, bbox.x0, bbox.y0, bbox.y1, bbox.z0, bbox.z1);
    Box3D domain2(bbox.x1, bbox.x1, bbox.y0, bbox.y1, bbox.z0, bbox.z1);
    Box3D domain3(bbox.x0 + 1, bbox.x1 - 1, bbox.y0, bbox.y0, bbox.z0, bbox.z1);
    Box3D domain4(bbox.x0 + 1, bbox.x1 - 1, bbox.y1, bbox.y1, bbox.z0, bbox.z1);
    Box3D domain5(bbox.x0 + 1, bbox.x1 - 1, bbox.y0 + 1, bbox.y1 - 1, bbox.z0, bbox.z0);
    Box3D domain6(bbox.x0 + 1, bbox.x1 - 1, bbox.y0 + 1, bbox.y1 - 1, bbox.z1, bbox.z1);
    Box3D intersection;
    NonLocalBoundaryDynamics3D<T, Descriptor> *dynamicsTemplate = 0;
    Dynamics<T, Descriptor> *baseDynamics = new NoDynamics<T, Descriptor>();
    switch (bcType) {
    case boundary::dirichlet:
        dynamicsTemplate = new NLD_VelocityBoundaryDynamics3D<T, Descriptor>(baseDynamics);
        break;
    case boundary::neumann:
        dynamicsTemplate = new NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>(baseDynamics);
        break;
    case boundary::freeslip: {
        bool noPenetration = true;
        dynamicsTemplate =
            new NLD_VelocityNeumannBoundaryDynamics3D<T, Descriptor>(baseDynamics, noPenetration);
    } break;
    case boundary::outflow:
        dynamicsTemplate = new NLD_OutflowBoundaryDynamics3D<T, Descriptor>(baseDynamics);
        break;
    default:
        PLB_ASSERT(0);  // Boundary type not implemented.
    }
    if (intersect(domain, domain1, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    if (intersect(domain, domain2, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    if (intersect(domain, domain3, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    if (intersect(domain, domain4, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    if (intersect(domain, domain5, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    if (intersect(domain, domain6, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    delete dynamicsTemplate;
}

template <typename T, template <typename U> class Descriptor>
void setAD_NLDboundaryDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, boundary::BcType bcType)
{
    setAD_NLDboundaryDynamics(lattice, lattice.getBoundingBox(), domain, bcType);
}

template <typename T, template <typename U> class Descriptor>
void setAD_NLDboundaryDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D bbox, Box3D domain, boundary::BcType bcType)
{
    Box3D domain1(bbox.x0, bbox.x0, bbox.y0, bbox.y1, bbox.z0, bbox.z1);
    Box3D domain2(bbox.x1, bbox.x1, bbox.y0, bbox.y1, bbox.z0, bbox.z1);
    Box3D domain3(bbox.x0 + 1, bbox.x1 - 1, bbox.y0, bbox.y0, bbox.z0, bbox.z1);
    Box3D domain4(bbox.x0 + 1, bbox.x1 - 1, bbox.y1, bbox.y1, bbox.z0, bbox.z1);
    Box3D domain5(bbox.x0 + 1, bbox.x1 - 1, bbox.y0 + 1, bbox.y1 - 1, bbox.z0, bbox.z0);
    Box3D domain6(bbox.x0 + 1, bbox.x1 - 1, bbox.y0 + 1, bbox.y1 - 1, bbox.z1, bbox.z1);
    Box3D intersection;
    NonLocalBoundaryDynamics3D<T, Descriptor> *dynamicsTemplate = 0;
    Dynamics<T, Descriptor> *baseDynamics = new NoDynamics<T, Descriptor>();
    switch (bcType) {
    case boundary::dirichlet:
        dynamicsTemplate = new NLD_AD_DirichletDynamics3D<T, Descriptor>(baseDynamics);
        break;
    case boundary::neumann:
        dynamicsTemplate = new NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>(baseDynamics);
        break;
    case boundary::freeslip:
        dynamicsTemplate = new NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>(baseDynamics);
        break;
    case boundary::outflow:
        dynamicsTemplate = new NLD_AD_NeumannBoundaryDynamics3D<T, Descriptor>(baseDynamics);
        break;
    default:
        PLB_ASSERT(0);  // Boundary type not implemented.
    }
    if (intersect(domain, domain1, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    if (intersect(domain, domain2, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    if (intersect(domain, domain3, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    if (intersect(domain, domain4, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    if (intersect(domain, domain5, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    if (intersect(domain, domain6, intersection)) {
        setNLDdynamics(lattice, intersection, dynamicsTemplate->clone());
    }
    delete dynamicsTemplate;
}

template <typename T, template <typename U> class Descriptor>
void instantiateOuterNLDboundary(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    instantiateOuterNLDboundary(lattice, lattice.getBoundingBox());
}

template <typename T, template <typename U> class Descriptor>
void instantiateOuterNLDboundary(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D bbox)
{
    BlockSurface3D surf(bbox, 1);
    integrateProcessingFunctional(
        new ExecutePlaneNLD_3D<T, Descriptor>(0, -1), surf.surface0N(), lattice);
    integrateProcessingFunctional(
        new ExecutePlaneNLD_3D<T, Descriptor>(0, +1), surf.surface0P(), lattice);
    integrateProcessingFunctional(
        new ExecutePlaneNLD_3D<T, Descriptor>(1, -1), surf.surface1N(), lattice);
    integrateProcessingFunctional(
        new ExecutePlaneNLD_3D<T, Descriptor>(1, +1), surf.surface1P(), lattice);
    integrateProcessingFunctional(
        new ExecutePlaneNLD_3D<T, Descriptor>(2, -1), surf.surface2N(), lattice);
    integrateProcessingFunctional(
        new ExecutePlaneNLD_3D<T, Descriptor>(2, +1), surf.surface2P(), lattice);

    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(0, -1, -1), surf.edge0NN(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(0, -1, 1), surf.edge0NP(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(0, 1, -1), surf.edge0PN(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(0, 1, 1), surf.edge0PP(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(1, -1, -1), surf.edge1NN(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(1, -1, 1), surf.edge1NP(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(1, 1, -1), surf.edge1PN(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(1, 1, 1), surf.edge1PP(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(2, -1, -1), surf.edge2NN(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(2, -1, 1), surf.edge2NP(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(2, 1, -1), surf.edge2PN(), lattice);
    integrateProcessingFunctional(
        new ExecuteEdgeNLD_3D<T, Descriptor>(2, 1, 1), surf.edge2PP(), lattice);

    integrateProcessingFunctional(
        new ExecuteCornerNLD_3D<T, Descriptor>(-1, -1, -1), surf.cornerNNN(), lattice);
    integrateProcessingFunctional(
        new ExecuteCornerNLD_3D<T, Descriptor>(-1, -1, 1), surf.cornerNNP(), lattice);
    integrateProcessingFunctional(
        new ExecuteCornerNLD_3D<T, Descriptor>(-1, 1, -1), surf.cornerNPN(), lattice);
    integrateProcessingFunctional(
        new ExecuteCornerNLD_3D<T, Descriptor>(-1, 1, 1), surf.cornerNPP(), lattice);
    integrateProcessingFunctional(
        new ExecuteCornerNLD_3D<T, Descriptor>(1, -1, -1), surf.cornerPNN(), lattice);
    integrateProcessingFunctional(
        new ExecuteCornerNLD_3D<T, Descriptor>(1, -1, 1), surf.cornerPNP(), lattice);
    integrateProcessingFunctional(
        new ExecuteCornerNLD_3D<T, Descriptor>(1, 1, -1), surf.cornerPPN(), lattice);
    integrateProcessingFunctional(
        new ExecuteCornerNLD_3D<T, Descriptor>(1, 1, 1), surf.cornerPPP(), lattice);
}

template <typename T, template <typename U> class Descriptor>
void instantiateOuterNLDboundary(Actions3D &action, plint blockNum, Box3D bbox)
{
    BlockSurface3D surf(bbox, 1);
    action.addProcessor(new ExecutePlaneNLD_3D<T, Descriptor>(0, -1), blockNum, surf.surface0N());
    action.addProcessor(new ExecutePlaneNLD_3D<T, Descriptor>(0, +1), blockNum, surf.surface0P());
    action.addProcessor(new ExecutePlaneNLD_3D<T, Descriptor>(1, -1), blockNum, surf.surface1N());
    action.addProcessor(new ExecutePlaneNLD_3D<T, Descriptor>(1, +1), blockNum, surf.surface1P());
    action.addProcessor(new ExecutePlaneNLD_3D<T, Descriptor>(2, -1), blockNum, surf.surface2N());
    action.addProcessor(new ExecutePlaneNLD_3D<T, Descriptor>(2, +1), blockNum, surf.surface2P());

    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(0, -1, -1), blockNum, surf.edge0NN());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(0, -1, 1), blockNum, surf.edge0NP());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(0, 1, -1), blockNum, surf.edge0PN());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(0, 1, 1), blockNum, surf.edge0PP());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(1, -1, -1), blockNum, surf.edge1NN());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(1, -1, 1), blockNum, surf.edge1NP());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(1, 1, -1), blockNum, surf.edge1PN());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(1, 1, 1), blockNum, surf.edge1PP());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(2, -1, -1), blockNum, surf.edge2NN());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(2, -1, 1), blockNum, surf.edge2NP());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(2, 1, -1), blockNum, surf.edge2PN());
    action.addProcessor(new ExecuteEdgeNLD_3D<T, Descriptor>(2, 1, 1), blockNum, surf.edge2PP());

    action.addProcessor(
        new ExecuteCornerNLD_3D<T, Descriptor>(-1, -1, -1), blockNum, surf.cornerNNN());
    action.addProcessor(
        new ExecuteCornerNLD_3D<T, Descriptor>(-1, -1, 1), blockNum, surf.cornerNNP());
    action.addProcessor(
        new ExecuteCornerNLD_3D<T, Descriptor>(-1, 1, -1), blockNum, surf.cornerNPN());
    action.addProcessor(
        new ExecuteCornerNLD_3D<T, Descriptor>(-1, 1, 1), blockNum, surf.cornerNPP());
    action.addProcessor(
        new ExecuteCornerNLD_3D<T, Descriptor>(1, -1, -1), blockNum, surf.cornerPNN());
    action.addProcessor(
        new ExecuteCornerNLD_3D<T, Descriptor>(1, -1, 1), blockNum, surf.cornerPNP());
    action.addProcessor(
        new ExecuteCornerNLD_3D<T, Descriptor>(1, 1, -1), blockNum, surf.cornerPPN());
    action.addProcessor(
        new ExecuteCornerNLD_3D<T, Descriptor>(1, 1, 1), blockNum, surf.cornerPPP());
}

template <typename T, template <class U> class Descriptor>
void setNLDdynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics)
{
    applyProcessingFunctional(
        new InstantiateNLDdynamicsFunctional3D<T, Descriptor>(dynamics), domain, lattice);
}

/* *************** Class InstantiateNLDdynamicsFunctional3D ************* */

template <typename T, template <typename U> class Descriptor>
InstantiateNLDdynamicsFunctional3D<T, Descriptor>::InstantiateNLDdynamicsFunctional3D(
    NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics_) :
    dynamics(dynamics_)
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateNLDdynamicsFunctional3D<T, Descriptor>::InstantiateNLDdynamicsFunctional3D(
    InstantiateNLDdynamicsFunctional3D<T, Descriptor> const &rhs) :
    dynamics(rhs.dynamics->clone())
{ }

template <typename T, template <typename U> class Descriptor>
InstantiateNLDdynamicsFunctional3D<T, Descriptor>
    &InstantiateNLDdynamicsFunctional3D<T, Descriptor>::operator=(
        InstantiateNLDdynamicsFunctional3D<T, Descriptor> const &rhs)
{
    delete dynamics;
    dynamics = rhs.dynamics->clone();
    return *this;
}

template <typename T, template <typename U> class Descriptor>
InstantiateNLDdynamicsFunctional3D<T, Descriptor>::~InstantiateNLDdynamicsFunctional3D()
{
    delete dynamics;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateNLDdynamicsFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                lattice.attributeDynamics(
                    iX, iY, iZ,
                    cloneAndInsertNLDatTop(
                        lattice.get(iX, iY, iZ).getDynamics(), dynamics->clone()));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT InstantiateNLDdynamicsFunctional3D<T, Descriptor>::appliesTo() const
{
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template <typename T, template <typename U> class Descriptor>
void InstantiateNLDdynamicsFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
InstantiateNLDdynamicsFunctional3D<T, Descriptor>
    *InstantiateNLDdynamicsFunctional3D<T, Descriptor>::clone() const
{
    return new InstantiateNLDdynamicsFunctional3D<T, Descriptor>(*this);
}

}  // namespace plb

#endif  // NLD_BOUNDARIES_3D_HH
