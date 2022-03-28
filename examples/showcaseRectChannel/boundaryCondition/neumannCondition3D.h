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
 * Neumann and outflow boundary conditions -- header file.
 */
#ifndef NEUMANN_CONDITION_3D_H
#define NEUMANN_CONDITION_3D_H

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "boundaryCondition/boundaryDynamics.h"
#include "core/blockLatticeBase3D.h"
#include "core/globalDefs.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class CopyUnknownPopulationsFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual CopyUnknownPopulationsFunctional3D<T, Descriptor, direction, orientation> *clone()
        const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
class CopyAllPopulationsFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual CopyAllPopulationsFunctional3D<T, Descriptor, normalX, normalY, normalZ> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
class CopyVelocityFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual CopyVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
class CopyTangentialVelocityFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual CopyTangentialVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ> *clone()
        const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
class CopyNormalVelocityFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual CopyNormalVelocityFunctional3D<T, Descriptor, normalX, normalY, normalZ> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor, int normalX, int normalY, int normalZ>
class CopyDensityFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual CopyDensityFunctional3D<T, Descriptor, normalX, normalY, normalZ> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/// Copy Unknown Populations Plus Impose Constant Pressure
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class FluidPressureOutlet3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual FluidPressureOutlet3D<T, Descriptor, direction, orientation> *clone() const
    {
        return new FluidPressureOutlet3D<T, Descriptor, direction, orientation>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
};

/// Copy Unknown Populations Plus Impose Constant Pressure
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class CopyUnknown3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual CopyUnknown3D<T, Descriptor, direction, orientation> *clone() const
    {
        return new CopyUnknown3D<T, Descriptor, direction, orientation>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }
};

/// Outflow boundary condition that works with external rhoBar and J from the previous time step.
///     WARNING 1: The provided rhoBar and J must contain the values of the
///     fields at the previous time step. This means that if  the
///      "external-rhoBar-J-collide-and-stream" is used at a processor level 0, then
///     this data processor must be integrated at processor level 1, and the data processor
///     to recompute the new rhoBar and J should be integrated at processor level 2.
///     WARNING 2: The lattice must be non-periodic in the direction of the outflow
///     normal. The reason for this is that the data processor uses data in the
///     envelope to access post-collision populations from the previous time step. if
///     the direction is periodic, these values get overwritten.
template <typename T, template <typename U> class Descriptor>
class VirtualOutlet : public BoxProcessingFunctional3D {
public:
    /* Type 0: Close to FluidPressureOutlet3D (imposes a strict pressure).
     * Type 1: Laplacian filter / extrapolation on the pressure / Free-surface like completion of
     *populations. Type 2: Laplacian filter / extrapolation on the pressure / Extrapolation of
     *non-equilibrium populations.
     **/
    VirtualOutlet(T outsideDensity_, Box3D globalDomain_, int type_ = 1);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual VirtualOutlet<T, Descriptor> *clone() const
    {
        return new VirtualOutlet<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;  // Block lattice.
        modified[1] = modif::nothing;          // RhoBar.
        modified[2] = modif::nothing;          // J.
    }

private:
    T outsideDensity;    // Boundary condition for the density (usually 1.0).
    Box3D globalDomain;  // The globalDomain must be at most as big as the whole simulation
                         // domain for non-periodic problems, and bigger than the whole simulation
                         // domain plus the envelope (per periodic direction) for periodic problems.
    int type;            // If type = 0 then this is very close to FluidPressureOutlet3D.
                         // If type = 1 some times gives the best results.
                         // If type = 2 is new and being tested.
};

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class VirtualOutletDynamics : public BoundaryCompositeDynamics<T, Descriptor> {
public:
    /// Constructor
    VirtualOutletDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision = true);
    VirtualOutletDynamics(HierarchicUnserializer &unserializer);
    /// Clone the object on its dynamic type.
    virtual VirtualOutletDynamics<T, Descriptor, direction, orientation> *clone() const;
    virtual void serialize(HierarchicSerializer &serializer) const;
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;
    void saveData(Cell<T, Descriptor> &cell) const;

private:
    mutable std::vector<T> savedFneq;
    mutable Array<T, Descriptor<T>::d> savedJ;
    mutable T savedRhoBar;
    static int id;
};

}  // namespace plb

#endif  // NEUMANN_CONDITION_3D_H
