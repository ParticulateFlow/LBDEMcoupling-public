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
 * local dynamics objects -- header file.
 */
#ifndef NLD_BOUNDARIES_3D_H
#define NLD_BOUNDARIES_3D_H

#include "boundaryCondition/NLD_boundaryDynamics3D.h"
#include "boundaryCondition/boundaryCondition.h"
#include "core/globalDefs.h"
#include "core/nonLocalDynamics3D.h"
#include "multiBlock/coupling3D.h"
#include "multiBlock/multiBlockLattice3D.h"

namespace plb {

/// A generic interface for non-local data processors that invoke the NLD dynamics objects.
/** These data processors don't do anything sophisticated. They simply call the dynamics
 *  objects with the right parameters, including the non-local information on neighboring
 *  cells.
 **/
template <typename T, template <typename U> class Descriptor>
class ExecuteNonLocalDynamics3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual ExecuteNonLocalDynamics3D<T, Descriptor> *clone() const;
    virtual int getStaticId() const
    {
        return staticId;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }

private:
    static const int staticId;
};

/// Invoke NLD dynamics objects on a plane surface.
template <typename T, template <typename U> class Descriptor>
class ExecutePlaneNLD_3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    ExecutePlaneNLD_3D();
    ExecutePlaneNLD_3D(int direction_, int orientation_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual ExecutePlaneNLD_3D<T, Descriptor> *clone() const;
    virtual int getStaticId() const
    {
        return staticId;
    }
    virtual void serialize(std::string &data) const;
    virtual void unserialize(std::string &data);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }

private:
    int direction, orientation;
    static const int staticId;
};

/// Invoke NLD dynamics objects on an edge.
template <typename T, template <typename U> class Descriptor>
class ExecuteEdgeNLD_3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    ExecuteEdgeNLD_3D();
    ExecuteEdgeNLD_3D(int plane_, int normal1_, int normal2_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual ExecuteEdgeNLD_3D<T, Descriptor> *clone() const;
    virtual int getStaticId() const
    {
        return staticId;
    }
    virtual void serialize(std::string &data) const;
    virtual void unserialize(std::string &data);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }

private:
    int plane, normal1, normal2;
    static const int staticId;
};

/// Invoke NLD dynamics objects on a corner.
template <typename T, template <typename U> class Descriptor>
class ExecuteCornerNLD_3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    ExecuteCornerNLD_3D();
    ExecuteCornerNLD_3D(int xNormal_, int yNormal_, int zNormal_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual ExecuteCornerNLD_3D<T, Descriptor> *clone() const;
    virtual int getStaticId() const
    {
        return staticId;
    }
    virtual void serialize(std::string &data) const;
    virtual void unserialize(std::string &data);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;
    }

private:
    int xNormal, yNormal, zNormal;
    static const int staticId;
};

/// Add the NLD data processors to the outer bounding-box of a lattice.
template <typename T, template <typename U> class Descriptor>
void instantiateOuterNLDboundary(MultiBlockLattice3D<T, Descriptor> &lattice);

/// Add the NLD data processors to the bounding-box of a given domain.
template <typename T, template <typename U> class Descriptor>
void instantiateOuterNLDboundary(MultiBlockLattice3D<T, Descriptor> &lattice, Box3D bbox);

/// Add the NLD data processors to the bounding-box of a given domain.
template <typename T, template <typename U> class Descriptor>
void instantiateOuterNLDboundary(Actions3D &action, plint blockNum, Box3D bbox);

/// Define the appropriate NLD dynamics, depending on whether it is
/// dirichlet, neumann, outflow, etc., on selected nodes of the bounding-box
/// of a given domain.
template <typename T, template <typename U> class Descriptor>
void setFluidNLDboundaryDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D bbox, Box3D domain, boundary::BcType bcType);

/// Define the appropriate NLD dynamics, depending on whether it is
/// dirichlet, neumann, outflow, etc., on the bounding-box of a given domain.
template <typename T, template <typename U> class Descriptor>
void setFluidNLDboundaryDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, boundary::BcType bcType);

/// Define the appropriate NLD dynamics, depending on whether it is
/// dirichlet, neumann, outflow, etc., on selected nodes of the bounding-box
/// of a given domain.
template <typename T, template <typename U> class Descriptor>
void setAD_NLDboundaryDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D bbox, Box3D domain, boundary::BcType bcType);

/// Define the appropriate NLD dynamics, depending on whether it is
/// dirichlet, neumann, outflow, etc., on the bounding-box of a given domain.
template <typename T, template <typename U> class Descriptor>
void setAD_NLDboundaryDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain, boundary::BcType bcType);

template <typename T, template <class U> class Descriptor>
void setNLDdynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, Box3D domain,
    NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics);

/* *************** Class InstantiateNLDdynamicsFunctional3D ************* */

template <typename T, template <typename U> class Descriptor>
class InstantiateNLDdynamicsFunctional3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    InstantiateNLDdynamicsFunctional3D(NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics_);
    InstantiateNLDdynamicsFunctional3D(
        InstantiateNLDdynamicsFunctional3D<T, Descriptor> const &rhs);
    InstantiateNLDdynamicsFunctional3D<T, Descriptor> &operator=(
        InstantiateNLDdynamicsFunctional3D<T, Descriptor> const &rhs);
    virtual ~InstantiateNLDdynamicsFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &lattice);
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual InstantiateNLDdynamicsFunctional3D<T, Descriptor> *clone() const;

private:
    NonLocalBoundaryDynamics3D<T, Descriptor> *dynamics;
};

}  // namespace plb

#endif  // NLD_BOUNDARIES_3D_H
