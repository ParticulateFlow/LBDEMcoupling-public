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

/** \file A helper for initialising 2D boundaries -- header file.  */

#ifndef BOUNDARY_CONDITION_2D_H
#define BOUNDARY_CONDITION_2D_H

#include "boundaryCondition/boundaryCondition.h"
#include "boundaryCondition/finiteDifferenceBoundaryProcessor2D.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class BlockLattice2D;
template <typename T, template <typename U> class Descriptor>
class MultiBlockLattice2D;

template <typename T, template <typename U> class Descriptor>
class OnLatticeBoundaryCondition2D {
public:
    virtual ~OnLatticeBoundaryCondition2D() { }
    virtual OnLatticeBoundaryCondition2D<T, Descriptor> *clone() const = 0;

    // PART I: Atomic-block version.

    virtual void addVelocityBoundary0N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addVelocityBoundary0P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addVelocityBoundary1N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addVelocityBoundary1P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    virtual void addPressureBoundary0N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addPressureBoundary0P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addPressureBoundary1N(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addPressureBoundary1P(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    virtual void addExternalVelocityCornerNN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addExternalVelocityCornerNP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addExternalVelocityCornerPN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addExternalVelocityCornerPP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    virtual void addInternalVelocityCornerNN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addInternalVelocityCornerNP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addInternalVelocityCornerPN(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addInternalVelocityCornerPP(
        plint x, plint y, BlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    /// Set velocity/Neumann condition on outer boundaries of the lattice (atomic-block
    ///   version).
    void setVelocityConditionOnBlockBoundaries(
        BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType = boundary::dirichlet);

    /// Set velocity/Neumann condition on a sub-domain, on the outer boundaries of
    ///   the lattice (atomic-block version).
    /** Attention: this function only has an effect when it is used on the outer surface
     *  of the atomic-block. For boundaries inside the domain, use
     *  the method which takes two Box2D arguments.
     **/
    void setVelocityConditionOnBlockBoundaries(
        BlockLattice2D<T, Descriptor> &lattice, Box2D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    /// Set velocity/Neumann condition on the block boundaries, but only on places which
    ///    intersect with the area of applicationDomain (atomic-block version).
    void setVelocityConditionOnBlockBoundaries(
        BlockLattice2D<T, Descriptor> &lattice, Box2D block, Box2D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    /// Set Pressure condition on outer boundaries of the lattice (atomic-block version).
    /** Attention: pressure conditions are implemented for edges only. On corners,
     *  this function has no effect.
     **/
    void setPressureConditionOnBlockBoundaries(
        BlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType = boundary::dirichlet);

    /// Set Pressure condition on a sub-domain, on the outer boundaries of
    ///   the lattice (atomic-block version).
    /** Attention: this function only has an effect when it is used on the outer surface
     *  of the atomic-block. For boundaries inside the domain, use
     *  the method which takes two Box2D arguments.
     *  Attention: pressure conditions are implemented for edges only. On corners,
     *  this function has no effect.
     **/
    void setPressureConditionOnBlockBoundaries(
        BlockLattice2D<T, Descriptor> &lattice, Box2D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    /// Set Pressure condition on the block boundaries, but only on places which
    ///    intersect with the area of applicationDomain (atomic-block version).
    /** Attention: pressure conditions are implemented for edges only. On corners,
     *  this function has no effect.
     **/
    void setPressureConditionOnBlockBoundaries(
        BlockLattice2D<T, Descriptor> &lattice, Box2D block, Box2D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    // PART II: Multi-block version.

    virtual void addVelocityBoundary0N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addVelocityBoundary0P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addVelocityBoundary1N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addVelocityBoundary1P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    virtual void addPressureBoundary0N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addPressureBoundary0P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addPressureBoundary1N(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addPressureBoundary1P(
        Box2D domain, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    virtual void addExternalVelocityCornerNN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addExternalVelocityCornerNP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addExternalVelocityCornerPN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addExternalVelocityCornerPP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    virtual void addInternalVelocityCornerNN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addInternalVelocityCornerNP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addInternalVelocityCornerPN(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;
    virtual void addInternalVelocityCornerPP(
        plint x, plint y, MultiBlockLattice2D<T, Descriptor> &lattice,
        boundary::BcType bcType = boundary::dirichlet) = 0;

    /// Set velocity/Neumann condition on outer boundaries of the lattice (multi-block
    ///   version).
    void setVelocityConditionOnBlockBoundaries(
        MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType = boundary::dirichlet);

    /// Set velocity/Neumann condition on a sub-domain, on the outer boundaries of
    ///   the lattice (multi-block version).
    /** Attention: this function only has an effect when it is used on the outer surface
     *  of the multi-block. For boundaries inside the domain, use
     *  the method which takes two Box2D arguments.
     **/
    void setVelocityConditionOnBlockBoundaries(
        MultiBlockLattice2D<T, Descriptor> &lattice, Box2D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    /// Set velocity/Neumann condition on the block boundaries, but only on places which
    ///    intersect with the area of applicationDomain (multi-block version).
    void setVelocityConditionOnBlockBoundaries(
        MultiBlockLattice2D<T, Descriptor> &lattice, Box2D block, Box2D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    /// Set Pressure condition on outer boundaries of the lattice (multi-block version).
    /** Attention: pressure conditions are implemented for edges only. On corners,
     *  this function has no effect.
     **/
    void setPressureConditionOnBlockBoundaries(
        MultiBlockLattice2D<T, Descriptor> &lattice, boundary::BcType bcType = boundary::dirichlet);

    /// Set Pressure condition on a sub-domain, on the outer boundaries of
    ///   the lattice (multi-block version).
    /** Attention: this function only has an effect when it is used on the outer surface
     *  of the multi-block. For boundaries inside the domain, use
     *  the method which takes two Box2D arguments.
     *  Attention: pressure conditions are implemented for edges only. On corners,
     *  this function has no effect.
     **/
    void setPressureConditionOnBlockBoundaries(
        MultiBlockLattice2D<T, Descriptor> &lattice, Box2D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    /// Set Pressure condition on the block boundaries, but only on places which
    ///    intersect with the area of applicationDomain (multi-block version).
    /** Attention: pressure conditions are implemented for edges only. On corners,
     *  this function has no effect.
     **/
    void setPressureConditionOnBlockBoundaries(
        MultiBlockLattice2D<T, Descriptor> &lattice, Box2D block, Box2D applicationDomain,
        boundary::BcType bcType = boundary::dirichlet);

    /// Set velocity/Neumann condition on inner boundaries of the lattice (multi-block
    ///   version).
    void setVelocityConditionOnInnerBlockBoundaries(
        MultiBlockLattice2D<T, Descriptor> &lattice, Box2D block,
        boundary::BcType bcType = boundary::dirichlet);
};

////////// Factory functions //////////////////////////////////////////////////

/// Generate the regularized boundary condition, managed by data processors.
template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createLocalBoundaryCondition2D();

/// Generate the regularized boundary condition, managed by dynamics only.
template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createDynamicsBasedLocalBoundaryCondition2D();

/// Generate a boundary condition which imposes equilibrium, but computes density properly.
template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createEquilibriumBoundaryCondition2D();

/// Generate Skordos boundary condition.
template <typename T, template <typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T, Descriptor> *createInterpBoundaryCondition2D();

}  // namespace plb

#endif  // BOUNDARY_CONDITION_2D_H
