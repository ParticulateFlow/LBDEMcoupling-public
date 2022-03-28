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
 * Helper functions for data field initialization -- generic implementation.
 */
#ifndef DATA_INITIALIZER_GENERICS_2D_H
#define DATA_INITIALIZER_GENERICS_2D_H

#include "atomicBlock/dataProcessorWrapper2D.h"
#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiGrid/multiGridUtil.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Initialization of the block-lattice *************** */
/* ******************************************************************* */

/* ************ Class SetCustomBoundaryVelocityFunctional2D ********** */

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
SetCustomBoundaryVelocityFunctional2D<
    T, Descriptor, VelocityFunction>::SetCustomBoundaryVelocityFunctional2D(VelocityFunction f_) :
    f(f_), velocityScale((T)1)
{ }

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
SetCustomBoundaryVelocityFunctional2D<T, Descriptor, VelocityFunction>
    *SetCustomBoundaryVelocityFunctional2D<T, Descriptor, VelocityFunction>::clone() const
{
    return new SetCustomBoundaryVelocityFunctional2D<T, Descriptor, VelocityFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
void SetCustomBoundaryVelocityFunctional2D<T, Descriptor, VelocityFunction>::execute(
    plint iX, plint iY, Cell<T, Descriptor> &cell) const
{
    Array<T, Descriptor<T>::d> u;
    f(iX, iY, u);
    u[0] *= velocityScale;
    u[1] *= velocityScale;
    cell.defineVelocity(u);
}

template <typename T, template <typename U> class Descriptor, class VelocityFunction>
void SetCustomBoundaryVelocityFunctional2D<T, Descriptor, VelocityFunction>::setscale(
    int dxScale, int dtScale)
{
    int dimDx = 1;
    int dimDt = -1;
    velocityScale = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
}

/* ************ Class SetCustomBoundaryDensityFunctional2D ********** */

template <typename T, template <typename U> class Descriptor, class DensityFunction>
SetCustomBoundaryDensityFunctional2D<
    T, Descriptor, DensityFunction>::SetCustomBoundaryDensityFunctional2D(DensityFunction f_) :
    f(f_)
{ }

template <typename T, template <typename U> class Descriptor, class DensityFunction>
void SetCustomBoundaryDensityFunctional2D<T, Descriptor, DensityFunction>::execute(
    plint iX, plint iY, Cell<T, Descriptor> &cell) const
{
    // No rescaling needed: rho is scale invariant.
    T rho = f(iX, iY);
    cell.defineDensity(rho);
}

template <typename T, template <typename U> class Descriptor, class DensityFunction>
SetCustomBoundaryDensityFunctional2D<T, Descriptor, DensityFunction>
    *SetCustomBoundaryDensityFunctional2D<T, Descriptor, DensityFunction>::clone() const
{
    return new SetCustomBoundaryDensityFunctional2D<T, Descriptor, DensityFunction>(*this);
}

/* ************ Class SetCustomBoundaryTemperatureFunctional2D ********** */

template <typename T, template <typename U> class Descriptor, class TemperatureFunction>
SetCustomBoundaryTemperatureFunctional2D<T, Descriptor, TemperatureFunction>::
    SetCustomBoundaryTemperatureFunctional2D(TemperatureFunction f_) :
    f(f_)
{ }

template <typename T, template <typename U> class Descriptor, class TemperatureFunction>
void SetCustomBoundaryTemperatureFunctional2D<T, Descriptor, TemperatureFunction>::execute(
    plint iX, plint iY, Cell<T, Descriptor> &cell) const
{
    // No rescaling needed: rho is scale invariant.
    T temperature = f(iX, iY);
    cell.defineTemperature(temperature);
}

template <typename T, template <typename U> class Descriptor, class TemperatureFunction>
SetCustomBoundaryTemperatureFunctional2D<T, Descriptor, TemperatureFunction>
    *SetCustomBoundaryTemperatureFunctional2D<T, Descriptor, TemperatureFunction>::clone() const
{
    return new SetCustomBoundaryTemperatureFunctional2D<T, Descriptor, TemperatureFunction>(*this);
}

/* ************ Class IniCustomEquilibriumFunctional2D ********** */

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
IniCustomEquilibriumFunctional2D<T, Descriptor, RhoUFunction>::IniCustomEquilibriumFunctional2D(
    RhoUFunction f_) :
    f(f_), velocityScale((T)1)
{ }

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomEquilibriumFunctional2D<T, Descriptor, RhoUFunction>::execute(
    plint iX, plint iY, Cell<T, Descriptor> &cell) const
{
    Array<T, Descriptor<T>::d> j;
    T rho;
    f(iX, iY, rho, j);
    Array<T, Descriptor<T>::d> force;
    force[0] = getExternalForceComponent(cell, 0);
    force[1] = getExternalForceComponent(cell, 1);
    j[0] = velocityScale * rho * (j[0] - (T)0.5 * force[0]);
    j[1] = velocityScale * rho * (j[1] - (T)0.5 * force[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T rhoBar = Descriptor<T>::rhoBar(rho);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
IniCustomEquilibriumFunctional2D<T, Descriptor, RhoUFunction>
    *IniCustomEquilibriumFunctional2D<T, Descriptor, RhoUFunction>::clone() const
{
    return new IniCustomEquilibriumFunctional2D<T, Descriptor, RhoUFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomEquilibriumFunctional2D<T, Descriptor, RhoUFunction>::setscale(
    int dxScale, int dtScale)
{
    int dimDx = 1;
    int dimDt = -1;
    velocityScale = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
}

/* ************ Class IniCustomThermalEquilibriumFunctional2D ********** */

template <typename T, template <typename U> class Descriptor, class RhoVelTempFunction>
IniCustomThermalEquilibriumFunctional2D<T, Descriptor, RhoVelTempFunction>::
    IniCustomThermalEquilibriumFunctional2D(RhoVelTempFunction f_) :
    f(f_), velocityScale((T)1)
{ }

template <typename T, template <typename U> class Descriptor, class RhoVelTempFunction>
void IniCustomThermalEquilibriumFunctional2D<T, Descriptor, RhoVelTempFunction>::execute(
    plint iX, plint iY, Cell<T, Descriptor> &cell) const
{
    Array<T, Descriptor<T>::d> j;
    T rho, temperature;
    f(iX, iY, rho, j, temperature);
    Array<T, Descriptor<T>::d> force;
    force[0] = getExternalForceComponent(cell, 0);
    force[1] = getExternalForceComponent(cell, 1);
    j[0] = velocityScale * rho * (j[0] - (T)0.5 * force[0]);
    j[1] = velocityScale * rho * (j[1] - (T)0.5 * force[1]);
    T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
    T thetaBar = temperature - (T)1;
    T rhoBar = Descriptor<T>::rhoBar(rho);
    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
        cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
    }
}

template <typename T, template <typename U> class Descriptor, class RhoVelTempFunction>
IniCustomThermalEquilibriumFunctional2D<T, Descriptor, RhoVelTempFunction>
    *IniCustomThermalEquilibriumFunctional2D<T, Descriptor, RhoVelTempFunction>::clone() const
{
    return new IniCustomThermalEquilibriumFunctional2D<T, Descriptor, RhoVelTempFunction>(*this);
}

template <typename T, template <typename U> class Descriptor, class RhoUFunction>
void IniCustomThermalEquilibriumFunctional2D<T, Descriptor, RhoUFunction>::setscale(
    int dxScale, int dtScale)
{
    int dimDx = 1;
    int dimDt = -1;
    velocityScale = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
}

/* ************ Class SetCustomOmegaFunctional2D ********** */

template <typename T, template <typename U> class Descriptor, class OmegaFunction>
SetCustomOmegaFunctional2D<T, Descriptor, OmegaFunction>::SetCustomOmegaFunctional2D(
    OmegaFunction f_) :
    f(f_)
{ }

template <typename T, template <typename U> class Descriptor, class OmegaFunction>
void SetCustomOmegaFunctional2D<T, Descriptor, OmegaFunction>::execute(
    plint iX, plint iY, Cell<T, Descriptor> &cell) const
{
    T omega = f(iX, iY);
    cell.getDynamics().setOmega(omega);
}

template <typename T, template <typename U> class Descriptor, class OmegaFunction>
SetCustomOmegaFunctional2D<T, Descriptor, OmegaFunction>
    *SetCustomOmegaFunctional2D<T, Descriptor, OmegaFunction>::clone() const
{
    return new SetCustomOmegaFunctional2D<T, Descriptor, OmegaFunction>(*this);
}

/* *************** PART II ******************************************* */
/* *************** Initialization of the scalar- and tensor-field **** */
/* ******************************************************************* */

/* ************ SetToScalarFunctionFunctional2D ********************** */

template <typename T, class Function>
SetToScalarFunctionFunctional2D<T, Function>::SetToScalarFunctionFunctional2D(Function f_) : f(f_)
{ }

template <typename T, class Function>
void SetToScalarFunctionFunctional2D<T, Function>::process(Box2D domain, ScalarField2D<T> &field)
{
    Dot2D relativeOffset = field.getLocation();
    Array<plint, 2> ofs(relativeOffset.x, relativeOffset.y);
    Array<plint, 2> pos;
    for (pos[0] = domain.x0; pos[0] <= domain.x1; ++pos[0]) {
        for (pos[1] = domain.y0; pos[1] <= domain.y1; ++pos[1]) {
            field.get(pos[0], pos[1]) = f(pos[0] + ofs[0], pos[1] + ofs[1]);
        }
    }
}

template <typename T, class Function>
SetToScalarFunctionFunctional2D<T, Function> *SetToScalarFunctionFunctional2D<T, Function>::clone()
    const
{
    return new SetToScalarFunctionFunctional2D<T, Function>(*this);
}

template <typename T, class Function>
BlockDomain::DomainT SetToScalarFunctionFunctional2D<T, Function>::appliesTo() const
{
    // Boundary cannot be included, because periodic boundaries would
    //   get the wrong value.
    return BlockDomain::bulk;
}

template <typename T, class Function>
void SetToScalarFunctionFunctional2D<T, Function>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************ SetToTensorFunctionFunctional2D ********************** */

template <typename T, int nDim, class Function>
SetToTensorFunctionFunctional2D<T, nDim, Function>::SetToTensorFunctionFunctional2D(Function f_) :
    f(f_)
{ }

template <typename T, int nDim, class Function>
void SetToTensorFunctionFunctional2D<T, nDim, Function>::process(
    Box2D domain, TensorField2D<T, nDim> &field)
{
    Dot2D relativeOffset = field.getLocation();
    Array<plint, 2> ofs(relativeOffset.x, relativeOffset.y);
    Array<plint, 2> pos;
    Array<T, nDim> value;
    for (pos[0] = domain.x0; pos[0] <= domain.x1; ++pos[0]) {
        for (pos[1] = domain.y0; pos[1] <= domain.y1; ++pos[1]) {
            f(pos[0] + ofs[0], pos[1] + ofs[1], value);
            field.get(pos[0], pos[1]) = value;
        }
    }
}

template <typename T, int nDim, class Function>
SetToTensorFunctionFunctional2D<T, nDim, Function>
    *SetToTensorFunctionFunctional2D<T, nDim, Function>::clone() const
{
    return new SetToTensorFunctionFunctional2D<T, nDim, Function>(*this);
}

template <typename T, int nDim, class Function>
BlockDomain::DomainT SetToTensorFunctionFunctional2D<T, nDim, Function>::appliesTo() const
{
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template <typename T, int nDim, class Function>
void SetToTensorFunctionFunctional2D<T, nDim, Function>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

}  // namespace plb

#endif  // DATA_INITIALIZER_GENERICS_2D_H
