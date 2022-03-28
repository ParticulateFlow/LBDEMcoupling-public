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

#ifndef ADVECTION_DIFFUSION_BOUNDARIES_HH
#define ADVECTION_DIFFUSION_BOUNDARIES_HH

#include "complexDynamics/advectionDiffusionBoundaries.h"
#include "complexDynamics/utilAdvectionDiffusion.h"
#include "core/util.h"
#include "latticeBoltzmann/advectionDiffusionDynamicsTemplates.h"
#include "latticeBoltzmann/advectionDiffusionLattices.h"
#include "latticeBoltzmann/indexTemplates.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void DensityClosure(Cell<T, Descriptor> &cell, Dynamics<T, Descriptor> const &dynamics)
{
    typedef Descriptor<T> D;

    T rho = dynamics.computeDensity(cell);
    T rhoBar = D::rhoBar(rho);

    plint missingNormal = 0;
    std::vector<plint> missingDiagonal =
        indexTemplates::subIndexOutgoing<D, direction, orientation>();
    std::vector<plint> knownIndexes = indexTemplates::remainingIndexes<D>(missingDiagonal);
    // here I know all missing and non missing f_i
    for (pluint iPop = 0; iPop < missingDiagonal.size(); ++iPop) {
        plint numOfNonNullComp = 0;
        for (int iDim = 0; iDim < D::d; ++iDim)
            numOfNonNullComp += abs(D::c[missingDiagonal[iPop]][iDim]);

        if (numOfNonNullComp == 1) {
            missingNormal = missingDiagonal[iPop];
            missingDiagonal.erase(missingDiagonal.begin() + iPop);
            break;
        }
    }

    T sum = T();
    for (pluint iPop = 0; iPop < knownIndexes.size(); ++iPop) {
        sum += cell[knownIndexes[iPop]];
    }
    cell[missingNormal] = rhoBar - sum;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void RegularizedClosure(Cell<T, Descriptor> &cell, Dynamics<T, Descriptor> const &dynamics)
{
    typedef Descriptor<T> D;
    typedef advectionDiffusionDynamicsTemplates<T, Descriptor> adTempl;

    T rhoBar = dynamics.computeRhoBar(cell);
    T *u = cell.getExternal(D::ExternalField::velocityBeginsAt);
    Array<T, D::d> jEq;
    for (plint iD = 0; iD < D::d; ++iD) {
        jEq[iD] = D::fullRho(rhoBar) * u[iD];
    }
    plint missingNormal = 0;
    std::vector<plint> missingDiagonal =
        indexTemplates::subIndexOutgoing<D, direction, orientation>();
    // here I know all missing and non missing f_i
    for (pluint iPop = 0; iPop < missingDiagonal.size(); ++iPop) {
        plint numOfNonNullComp = 0;
        for (int iDim = 0; iDim < D::d; ++iDim)
            numOfNonNullComp += abs(D::c[missingDiagonal[iPop]][iDim]);

        if (numOfNonNullComp == 1) {
            missingNormal = missingDiagonal[iPop];
            missingDiagonal.erase(missingDiagonal.begin() + iPop);
            break;
        }
    }

    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices

    cell[missingNormal] =
        adTempl::bgk_ma1_equilibrium(missingNormal, rhoBar, jEq)
        - (cell[indexTemplates::opposite<D>(missingNormal)]
           - adTempl::bgk_ma1_equilibrium(indexTemplates::opposite<D>(missingNormal), rhoBar, jEq));
}

// ============= flat wall standard boundary ==================//

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_AdvectionDiffusion") + util::val2str(direction) + std::string("_")
        + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::
    AdvectionDiffusionBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::
    AdvectionDiffusionBoundaryDynamics(HierarchicUnserializer &unserializer) :
    StoreDensityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>
    *AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::clone() const
{
    return new AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    DensityClosure<T, Descriptor, direction, orientation>(cell, *this);
}

// ============= flat wall regularized boundary ==================//

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor,
        RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation> >(
        std::string("Boundary_RegularizedAdvectionDiffusion") + util::val2str(direction)
        + std::string("_") + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::
    RegularizedAdvectionDiffusionBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::
    RegularizedAdvectionDiffusionBoundaryDynamics(HierarchicUnserializer &unserializer) :
    StoreDensityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>
    *RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::clone()
        const
{
    return new RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>(
        *this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::getId()
    const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void RegularizedAdvectionDiffusionBoundaryDynamics<
    T, Descriptor, direction, orientation>::completePopulations(Cell<T, Descriptor> &cell) const
{
    // RegularizedClosure<T,Descriptor,direction,orientation>(cell, *this);

    typedef Descriptor<T> D;

    T rhoBar = this->computeRhoBar(cell);

    Array<T, D::d> jEq, jNeq;
    jEq.from_cArray(cell.getExternal(D::ExternalField::velocityBeginsAt));
    jEq *= D::fullRho(rhoBar);

    T jEqSqr = VectorTemplate<T, Descriptor>::normSqr(jEq);
    boundaryTemplates<T, Descriptor, direction, orientation>::compute_jNeq(
        cell.getDynamics(), cell, rhoBar, jEq, jEqSqr, jNeq);

    Array<T, SymmetricTensor<T, Descriptor>::n> dummyPiNeq;
    cell.getDynamics().regularize(cell, rhoBar, jEq + jNeq, T(), dummyPiNeq);
}

// ============= flat wall complete regularized boundary ==================//

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int CompleteRegularizedAdvectionDiffusionBoundaryDynamics<
    T, Descriptor, direction, orientation>::id =
    meta::registerGeneralDynamics<
        T, Descriptor,
        CompleteRegularizedAdvectionDiffusionBoundaryDynamics<
            T, Descriptor, direction, orientation> >(
        std::string("Boundary_CompleteRegularizedAdvectionDiffusion") + util::val2str(direction)
        + std::string("_") + util::val2str(orientation));

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
CompleteRegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::
    CompleteRegularizedAdvectionDiffusionBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
CompleteRegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>::
    CompleteRegularizedAdvectionDiffusionBoundaryDynamics(HierarchicUnserializer &unserializer) :
    StoreDensityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
CompleteRegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>
    *CompleteRegularizedAdvectionDiffusionBoundaryDynamics<
        T, Descriptor, direction, orientation>::clone() const
{
    return new CompleteRegularizedAdvectionDiffusionBoundaryDynamics<
        T, Descriptor, direction, orientation>(*this);
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
int CompleteRegularizedAdvectionDiffusionBoundaryDynamics<
    T, Descriptor, direction, orientation>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int direction, int orientation>
void CompleteRegularizedAdvectionDiffusionBoundaryDynamics<
    T, Descriptor, direction, orientation>::completePopulations(Cell<T, Descriptor> &cell) const
{
    typedef Descriptor<T> D;

    T phiBar = this->computeRhoBar(cell);

    Array<T, D::d> jEq, jNeq;
    jEq.from_cArray(cell.getExternal(D::ExternalField::velocityBeginsAt));
    T rhoBar = *cell.getExternal(D::ExternalField::rhoBarBeginsAt);
    T rho = D::fullRho(rhoBar);
    T phi = D::fullRho(phiBar);

    jEq *= phi * rho;

    T rhoPhiBar = D::rhoBar(rho * phi);
    T jEqSqr = VectorTemplate<T, Descriptor>::normSqr(jEq);
    boundaryTemplates<T, Descriptor, direction, orientation>::compute_jNeq(
        cell.getDynamics(), cell, rhoPhiBar, jEq, jEqSqr, jNeq);

    Array<T, SymmetricTensor<T, Descriptor>::n> dummyPiNeq;
    cell.getDynamics().regularize(cell, rhoPhiBar, jEq + jNeq, T(), dummyPiNeq);
}

// =============== 2D corners ===================//

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
int AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal> >(
        std::string("Boundary_AdvectionDiffusionCorner") + util::val2str(xNormal) + std::string("_")
        + util::val2str(yNormal));

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal>::
    AdvectionDiffusionCornerDynamics2D(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal>::
    AdvectionDiffusionCornerDynamics2D(HierarchicUnserializer &unserializer) :
    StoreDensityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal>
    *AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal>::clone() const
{
    return new AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal>(*this);
}

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
int AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
void AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    typedef Descriptor<T> D;
    typedef advectionDiffusionDynamicsTemplates<T, Descriptor> adTempl;

    T rhoBar = this->computeRhoBar(cell);
    T *u = cell.getExternal(D::ExternalField::velocityBeginsAt);
    Array<T, D::d> jEq;
    for (plint iD = 0; iD < D::d; ++iD) {
        jEq[iD] = D::fullRho(rhoBar) * u[iD];
    }
    // I need to get Missing information on the corners !!!!
    std::vector<plint> unknownIndexes =
        utilAdvDiff::subIndexOutgoing2DonCorners<D, xNormal, yNormal>();
    // here I know all missing and non missing f_i

    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices

    for (pluint iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
        cell[unknownIndexes[iPop]] =
            adTempl::bgk_ma1_equilibrium(unknownIndexes[iPop], rhoBar, jEq)
            - (cell[indexTemplates::opposite<D>(unknownIndexes[iPop])]
               - adTempl::bgk_ma1_equilibrium(
                   indexTemplates::opposite<D>(unknownIndexes[iPop]), rhoBar, jEq));
    }
}

// =============== 3D corners ===================//

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
int AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::id =
    meta::registerGeneralDynamics<
        T, Descriptor,
        AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal> >(
        std::string("Boundary_AdvectionDiffusionCorner") + util::val2str(xNormal) + std::string("_")
        + util::val2str(yNormal) + std::string("_") + util::val2str(zNormal));

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::
    AdvectionDiffusionCornerDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::
    AdvectionDiffusionCornerDynamics3D(HierarchicUnserializer &unserializer) :
    StoreDensityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>
    *AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::clone() const
{
    return new AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>(*this);
}

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
int AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionCornerDynamics3D<
    T, Descriptor, xNormal, yNormal, zNormal>::completePopulations(Cell<T, Descriptor> &cell) const
{
    typedef Descriptor<T> D;
    typedef advectionDiffusionDynamicsTemplates<T, Descriptor> adTempl;

    T rhoBar = this->computeRhoBar(cell);
    T *u = cell.getExternal(D::ExternalField::velocityBeginsAt);
    Array<T, D::d> jEq;
    for (plint iD = 0; iD < D::d; ++iD) {
        jEq[iD] = D::fullRho(rhoBar) * u[iD];
    }
    // I need to get Missing information on the corners !!!!
    std::vector<plint> unknownIndexes =
        utilAdvDiff::subIndexOutgoing3DonCorners<D, xNormal, yNormal, zNormal>();
    // here I know all missing and non missing f_i

    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices

    for (pluint iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
        cell[unknownIndexes[iPop]] =
            adTempl::bgk_ma1_equilibrium(unknownIndexes[iPop], rhoBar, jEq)
            - (cell[indexTemplates::opposite<D>(unknownIndexes[iPop])]
               - adTempl::bgk_ma1_equilibrium(
                   indexTemplates::opposite<D>(unknownIndexes[iPop]), rhoBar, jEq));
    }
}

// =============== 3D edges ===================//

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
int AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::id =
    meta::registerGeneralDynamics<
        T, Descriptor, AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2> >(
        std::string("Boundary_AdvectionDiffusionEdge") + util::val2str(plane) + std::string("_")
        + util::val2str(normal1) + std::string("_") + util::val2str(normal2));

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::
    AdvectionDiffusionEdgeDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision) :
    StoreDensityDynamics<T, Descriptor>(baseDynamics, automaticPrepareCollision)
{ }

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::
    AdvectionDiffusionEdgeDynamics3D(HierarchicUnserializer &unserializer) :
    StoreDensityDynamics<T, Descriptor>(0, false)
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>
    *AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::clone() const
{
    return new AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>(*this);
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
int AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
void AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2>::completePopulations(
    Cell<T, Descriptor> &cell) const
{
    typedef Descriptor<T> D;
    typedef advectionDiffusionDynamicsTemplates<T, Descriptor> adTempl;

    T rhoBar = this->computeRhoBar(cell);
    T *u = cell.getExternal(D::ExternalField::velocityBeginsAt);
    Array<T, D::d> jEq;
    for (plint iD = 0; iD < D::d; ++iD) {
        jEq[iD] = D::fullRho(rhoBar) * u[iD];
    }
    // I need to get Missing information on the corners !!!!
    std::vector<plint> unknownIndexes =
        utilAdvDiff::subIndexOutgoing3DonEdges<D, plane, normal1, normal2>();
    // here I know all missing and non missing f_i

    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices

    for (pluint iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
        cell[unknownIndexes[iPop]] =
            adTempl::bgk_ma1_equilibrium(unknownIndexes[iPop], rhoBar, jEq)
            - (cell[indexTemplates::opposite<D>(unknownIndexes[iPop])]
               - adTempl::bgk_ma1_equilibrium(
                   indexTemplates::opposite<D>(unknownIndexes[iPop]), rhoBar, jEq));
    }
}

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_BOUNDARIES_HH
