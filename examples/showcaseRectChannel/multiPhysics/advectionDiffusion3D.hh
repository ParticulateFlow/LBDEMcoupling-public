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

#ifndef ADVECTION_DIFFUSION_3D_HH
#define ADVECTION_DIFFUSION_3D_HH

#include <cmath>

#include "atomicBlock/blockLattice3D.h"
#include "core/dynamics.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiPhysics/advectionDiffusion3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"

namespace plb {

template <typename T, template <typename U> class TemperatureDescriptor>
void VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>::process(
    Box3D domain, BlockLattice3D<T, TemperatureDescriptor> &temperature,
    TensorField3D<T, 3> &velocity)
{
    Dot3D offset = computeRelativeDisplacement(temperature, velocity);
    const int velOffset = TemperatureDescriptor<T>::ExternalField::velocityBeginsAt;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T *u = temperature.get(iX, iY, iZ).getExternal(velOffset);
                velocity.get(iX + offset.x, iY + offset.y, iZ + offset.z).to_cArray(u);
            }
        }
    }
}

template <typename T, template <typename U> class TemperatureDescriptor>
VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>
    *VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>::clone() const
{
    return new VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>(*this);
}

template <typename T, template <typename U> class TemperatureDescriptor>
BlockDomain::DomainT VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class TemperatureDescriptor>
void VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class TemperatureDescriptor>
void velocityToPassiveAdvDiff(
    MultiBlockLattice3D<T, TemperatureDescriptor> &temperature, MultiTensorField3D<T, 3> &velocity,
    Box3D domain)
{
    applyProcessingFunctional(
        new VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>(), domain, temperature, velocity);
}

template <typename T, template <typename U> class TemperatureDescriptor>
void N_VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>::process(
    Box3D domain, BlockLattice3D<T, TemperatureDescriptor> &temperature,
    NTensorField3D<T> &velocity)
{
    PLB_PRECONDITION(velocity.getNdim() == 3);
    Dot3D offset = computeRelativeDisplacement(temperature, velocity);
    const int velOffset = TemperatureDescriptor<T>::ExternalField::velocityBeginsAt;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T *u_to = temperature.get(iX, iY, iZ).getExternal(velOffset);
                T *u_from = velocity.get(iX + offset.x, iY + offset.y, iZ + offset.z);
                u_to[0] = u_from[0];
                u_to[1] = u_from[1];
                u_to[2] = u_from[2];
            }
        }
    }
}

template <typename T, template <typename U> class TemperatureDescriptor>
N_VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>
    *N_VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>::clone() const
{
    return new N_VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>(*this);
}

template <typename T, template <typename U> class TemperatureDescriptor>
BlockDomain::DomainT N_VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class TemperatureDescriptor>
void N_VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class TemperatureDescriptor>
void NVelocityToPassiveAdvDiff(
    MultiBlockLattice3D<T, TemperatureDescriptor> &temperature, MultiNTensorField3D<T> &velocity,
    Box3D domain)
{
    PLB_PRECONDITION(velocity.getNdim() == 3);
    applyProcessingFunctional(
        new N_VelocityToPassiveAdvDiff3D<T, TemperatureDescriptor>(), domain, temperature,
        velocity);
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
LatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::LatticeToPassiveAdvDiff3D(
    T scaling_) :
    scaling(scaling_)
{ }

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void LatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::process(
    Box3D domain, BlockLattice3D<T, FluidDescriptor> &fluid,
    BlockLattice3D<T, ScalarDescriptor> &scalar)
{
    Dot3D offset = computeRelativeDisplacement(fluid, scalar);
    const int velOffset = ScalarDescriptor<T>::ExternalField::velocityBeginsAt;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T *u =
                    scalar.get(iX + offset.x, iY + offset.y, iZ + offset.z).getExternal(velOffset);
                Array<T, 3> velocity;
                fluid.get(iX, iY, iZ).computeVelocity(velocity);
                velocity *= scaling;
                velocity.to_cArray(u);
            }
        }
    }
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
LatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>
    *LatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::clone() const
{
    return new LatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>(*this);
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
BlockDomain::DomainT LatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::appliesTo()
    const
{
    return BlockDomain::bulk;
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void LatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void latticeToPassiveAdvDiff(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, ScalarDescriptor> &scalar, Box3D domain)
{
    applyProcessingFunctional(
        new LatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>(), domain, fluid,
        scalar);
}

// =====  complete regularized advection diffusion coupling implementation
template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
LatticeToPassiveComplRegAdvDiff3D<
    T, FluidDescriptor, ScalarDescriptor>::LatticeToPassiveComplRegAdvDiff3D(T scaling_) :
    scaling(scaling_)
{ }

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void LatticeToPassiveComplRegAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::process(
    Box3D domain, BlockLattice3D<T, FluidDescriptor> &fluid,
    BlockLattice3D<T, ScalarDescriptor> &scalar)
{
    Dot3D offset = computeRelativeDisplacement(fluid, scalar);

    enum {
        rhoBarOffset = ScalarDescriptor<T>::ExternalField::rhoBarBeginsAt,
        velocityOffset = ScalarDescriptor<T>::ExternalField::velocityBeginsAt,
        piNeqOffset = ScalarDescriptor<T>::ExternalField::piNeqBeginsAt,
        omegaOffset = ScalarDescriptor<T>::ExternalField::omegaBeginsAt
    };

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // T *u = scalar.get(iX+offset.x,iY+offset.y,iZ+offset.z).getExternal(velOffset);
                // Array<T,3> velocity;
                // fluid.get(iX,iY,iZ).computeVelocity(velocity);
                // velocity *= scaling;
                // velocity.to_cArray(u);

                // Velocity coupling
                T rhoBar;
                Array<T, FluidDescriptor<T>::d> j;
                Array<T, SymmetricTensor<T, FluidDescriptor>::n> piNeq;
                fluid.get(iX, iY, iZ)
                    .getDynamics()
                    .computeRhoBarJPiNeq(fluid.get(iX, iY, iZ), rhoBar, j, piNeq);

                Array<T, FluidDescriptor<T>::d> vel =
                    FluidDescriptor<T>::invRho(rhoBar) * j * scaling;

                *scalar.get(iX + offset.x, iY + offset.y, iZ + offset.z).getExternal(omegaOffset) =
                    fluid.get(iX, iY, iZ).getDynamics().getOmega();
                *scalar.get(iX + offset.x, iY + offset.y, iZ + offset.z).getExternal(rhoBarOffset) =
                    rhoBar;
                vel.to_cArray(scalar.get(iX + offset.x, iY + offset.y, iZ + offset.z)
                                  .getExternal(velocityOffset));
                piNeq.to_cArray(scalar.get(iX + offset.x, iY + offset.y, iZ + offset.z)
                                    .getExternal(piNeqOffset));
            }
        }
    }
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
LatticeToPassiveComplRegAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>
    *LatticeToPassiveComplRegAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::clone() const
{
    return new LatticeToPassiveComplRegAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>(*this);
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
BlockDomain::DomainT
    LatticeToPassiveComplRegAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void LatticeToPassiveComplRegAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void latticeToPassiveComplRegAdvDiff(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, ScalarDescriptor> &scalar, Box3D domain)
{
    applyProcessingFunctional(
        new LatticeToPassiveComplRegAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>(), domain,
        fluid, scalar);
}

// Implementation of lattice to passive coupling for turbulent flows.

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
TurbulentLatticeToPassiveAdvDiff3D<
    T, FluidDescriptor, ScalarDescriptor>::TurbulentLatticeToPassiveAdvDiff3D(T Pr_t)
{
    PLB_ASSERT(!util::isZero(Pr_t));
    C = FluidDescriptor<T>::cs2 * ScalarDescriptor<T>::invCs2 / Pr_t;
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void TurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::process(
    Box3D domain, BlockLattice3D<T, FluidDescriptor> &fluid,
    BlockLattice3D<T, ScalarDescriptor> &scalar)
{
    Dot3D offset = computeRelativeDisplacement(fluid, scalar);
    const int velOffset = ScalarDescriptor<T>::ExternalField::velocityBeginsAt;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, FluidDescriptor> &fluidCell = fluid.get(iX, iY, iZ);
                Cell<T, ScalarDescriptor> &scalarCell =
                    scalar.get(iX + offset.x, iY + offset.y, iZ + offset.z);

                // Set the advection velocity.
                T *u = scalarCell.getExternal(velOffset);
                Array<T, 3> velocity;
                fluidCell.computeVelocity(velocity);
                velocity.to_cArray(u);

                // Set the relaxation parameter for the advection-diffusion.
                T fluidOmega = fluidCell.getDynamics().getDynamicParameter(
                    dynamicParams::dynamicOmega, fluidCell);
                if (!util::isZero(fluidOmega)) {
                    T fluidOmega0 = fluidCell.getDynamics().getOmega();
                    T scalarOmega0 = scalarCell.getDynamics().getOmega();
                    T scalarOmega =
                        (T)1 / ((T)1 / scalarOmega0 + C * ((T)1 / fluidOmega - (T)1 / fluidOmega0));
                    scalarCell.getDynamics().setOmega(scalarOmega);
                }
            }
        }
    }
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
TurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>
    *TurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::clone() const
{
    return new TurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>(*this);
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
BlockDomain::DomainT
    TurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void TurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;       // Fluid lattice.
    modified[1] = modif::allVariables;  // Advection-Diffusion lattice.
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void turbulentLatticeToPassiveAdvDiff(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, ScalarDescriptor> &scalar, T Pr_t, Box3D domain)
{
    applyProcessingFunctional(
        new TurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>(Pr_t), domain,
        fluid, scalar);
}

// Implementation of a masked lattice to passive coupling for turbulent flows.

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::
    MaskedTurbulentLatticeToPassiveAdvDiff3D(T Pr_t, int whichFlag_) :
    whichFlag(whichFlag_)
{
    PLB_ASSERT(!util::isZero(Pr_t));
    C = FluidDescriptor<T>::cs2 * ScalarDescriptor<T>::invCs2 / Pr_t;
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::
    processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    BlockLattice3D<T, FluidDescriptor> *fluid =
        dynamic_cast<BlockLattice3D<T, FluidDescriptor> *>(blocks[0]);
    BlockLattice3D<T, ScalarDescriptor> *scalar =
        dynamic_cast<BlockLattice3D<T, ScalarDescriptor> *>(blocks[1]);
    ScalarField3D<int> *mask = dynamic_cast<ScalarField3D<int> *>(blocks[2]);
    PLB_ASSERT(fluid);
    PLB_ASSERT(scalar);
    PLB_ASSERT(mask);

    Dot3D ofsScalar = computeRelativeDisplacement(*fluid, *scalar);
    Dot3D ofsMask = computeRelativeDisplacement(*fluid, *mask);
    const int velOffset = ScalarDescriptor<T>::ExternalField::velocityBeginsAt;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask->get(iX + ofsMask.x, iY + ofsMask.y, iZ + ofsMask.z) == whichFlag) {
                    Cell<T, FluidDescriptor> &fluidCell = fluid->get(iX, iY, iZ);
                    Cell<T, ScalarDescriptor> &scalarCell =
                        scalar->get(iX + ofsScalar.x, iY + ofsScalar.y, iZ + ofsScalar.z);

                    // Set the advection velocity.
                    T *u = scalarCell.getExternal(velOffset);
                    Array<T, 3> velocity;
                    fluidCell.computeVelocity(velocity);
                    velocity.to_cArray(u);

                    // Set the relaxation parameter for the advection-diffusion.
                    T fluidOmega = fluidCell.getDynamics().getDynamicParameter(
                        dynamicParams::dynamicOmega, fluidCell);
                    if (!util::isZero(fluidOmega)) {
                        T fluidOmega0 = fluidCell.getDynamics().getOmega();
                        T scalarOmega0 = scalarCell.getDynamics().getOmega();
                        T scalarOmega =
                            (T)1
                            / ((T)1 / scalarOmega0 + C * ((T)1 / fluidOmega - (T)1 / fluidOmega0));
                        scalarCell.getDynamics().setOmega(scalarOmega);
                    }
                }
            }
        }
    }
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>
    *MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::clone() const
{
    return new MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>(
        *this);
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
BlockDomain::DomainT MaskedTurbulentLatticeToPassiveAdvDiff3D<
    T, FluidDescriptor, ScalarDescriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;       // Fluid lattice.
    modified[1] = modif::allVariables;  // Advection-Diffusion lattice.
    modified[2] = modif::nothing;       // Mask.
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void maskedTurbulentLatticeToPassiveAdvDiff(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, ScalarDescriptor> &scalar, MultiScalarField3D<int> &mask, T Pr_t,
    int whichFlag, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&fluid);
    args.push_back(&scalar);
    args.push_back(&mask);
    applyProcessingFunctional(
        new MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>(
            Pr_t, whichFlag),
        domain, args);
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
N_MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::
    N_MaskedTurbulentLatticeToPassiveAdvDiff3D(T Pr_t, int whichFlag_) :
    whichFlag(whichFlag_)
{
    PLB_ASSERT(!util::isZero(Pr_t));
    C = FluidDescriptor<T>::cs2 * ScalarDescriptor<T>::invCs2 / Pr_t;
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void N_MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::
    processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    BlockLattice3D<T, FluidDescriptor> *fluid =
        dynamic_cast<BlockLattice3D<T, FluidDescriptor> *>(blocks[0]);
    BlockLattice3D<T, ScalarDescriptor> *scalar =
        dynamic_cast<BlockLattice3D<T, ScalarDescriptor> *>(blocks[1]);
    NTensorField3D<int> *mask = dynamic_cast<NTensorField3D<int> *>(blocks[2]);
    PLB_ASSERT(fluid);
    PLB_ASSERT(scalar);
    PLB_ASSERT(mask);

    Dot3D ofsScalar = computeRelativeDisplacement(*fluid, *scalar);
    Dot3D ofsMask = computeRelativeDisplacement(*fluid, *mask);
    const int velOffset = ScalarDescriptor<T>::ExternalField::velocityBeginsAt;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (*mask->get(iX + ofsMask.x, iY + ofsMask.y, iZ + ofsMask.z) == whichFlag) {
                    Cell<T, FluidDescriptor> &fluidCell = fluid->get(iX, iY, iZ);
                    Cell<T, ScalarDescriptor> &scalarCell =
                        scalar->get(iX + ofsScalar.x, iY + ofsScalar.y, iZ + ofsScalar.z);

                    // Set the advection velocity.
                    T *u = scalarCell.getExternal(velOffset);
                    Array<T, 3> velocity;
                    fluidCell.computeVelocity(velocity);
                    velocity.to_cArray(u);

                    // Set the relaxation parameter for the advection-diffusion.
                    T fluidOmega = fluidCell.getDynamics().getDynamicParameter(
                        dynamicParams::dynamicOmega, fluidCell);
                    if (!util::isZero(fluidOmega)) {
                        T fluidOmega0 = fluidCell.getDynamics().getOmega();
                        T scalarOmega0 = scalarCell.getDynamics().getOmega();
                        T scalarOmega =
                            (T)1
                            / ((T)1 / scalarOmega0 + C * ((T)1 / fluidOmega - (T)1 / fluidOmega0));
                        scalarCell.getDynamics().setOmega(scalarOmega);
                    }
                }
            }
        }
    }
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
N_MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>
    *N_MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::clone() const
{
    return new N_MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>(
        *this);
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
BlockDomain::DomainT N_MaskedTurbulentLatticeToPassiveAdvDiff3D<
    T, FluidDescriptor, ScalarDescriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void N_MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;       // Fluid lattice.
    modified[1] = modif::allVariables;  // Advection-Diffusion lattice.
    modified[2] = modif::nothing;       // Mask.
}

template <
    typename T, template <typename U1> class FluidDescriptor,
    template <typename U2> class ScalarDescriptor>
void NMaskedTurbulentLatticeToPassiveAdvDiff(
    MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiBlockLattice3D<T, ScalarDescriptor> &scalar, MultiNTensorField3D<int> &mask, T Pr_t,
    int whichFlag, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&fluid);
    args.push_back(&scalar);
    args.push_back(&mask);
    applyProcessingFunctional(
        new N_MaskedTurbulentLatticeToPassiveAdvDiff3D<T, FluidDescriptor, ScalarDescriptor>(
            Pr_t, whichFlag),
        domain, args);
}

template <typename T, template <typename U> class Descriptor>
CrystallizeAndAggregate<T, Descriptor>::CrystallizeAndAggregate(T Ncr_, T Nag_) :
    Ncr(Ncr_), Nag(Nag_)
{ }

template <typename T, template <typename U> class Descriptor>
void CrystallizeAndAggregate<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    plint bbId = BounceBack<T, Descriptor>().getId();
    Box3D extendedDomain(domain.enlarge(1));
    ScalarField3D<int> bbNodes(lattice.getNx(), lattice.getNy(), lattice.getNz());
    for (plint iX = extendedDomain.x0; iX <= extendedDomain.x1; ++iX) {
        for (plint iY = extendedDomain.y0; iY <= extendedDomain.y1; ++iY) {
            for (plint iZ = extendedDomain.z0; iZ <= extendedDomain.z1; ++iZ) {
                bbNodes.get(iX, iY, iZ) =
                    lattice.get(iX, iY, iZ).getDynamics().getId() == bbId ? 1 : 0;
            }
        }
    }
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (!bbNodes.get(iX, iY, iZ)) {
                    T N = lattice.get(iX, iY, iZ).computeDensity();
                    if (N >= Ncr) {
                        lattice.attributeDynamics(iX, iY, iZ, new BounceBack<T, Descriptor>());
                    } else if (N >= Nag) {
                        int numNeighbors = 0;
                        for (plint dx = -1; dx <= +1; ++dx) {
                            for (plint dy = -1; dy <= +1; ++dy) {
                                for (plint dz = -1; dz <= +1; ++dz) {
                                    numNeighbors += bbNodes.get(iX, iY, iZ);
                                }
                            }
                        }
                        if (numNeighbors > 0) {
                            lattice.attributeDynamics(iX, iY, iZ, new BounceBack<T, Descriptor>());
                        }
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CrystallizeAndAggregate<T, Descriptor> *CrystallizeAndAggregate<T, Descriptor>::clone() const
{
    return new CrystallizeAndAggregate<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT CrystallizeAndAggregate<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void CrystallizeAndAggregate<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dynamicVariables;
}

template <typename T, template <typename U> class Descriptor>
void crystallizeAndAggregate(
    MultiBlockLattice3D<T, Descriptor> &lattice, T Ncr, T Nag, Box3D domain)
{
    applyProcessingFunctional(
        new CrystallizeAndAggregate<T, Descriptor>(Ncr, Nag), domain, lattice);
}

/* Coupling between advection-diffusion and free-surface. */

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    AdvectionDiffusionFreeSurfaceCoupling3D(
        Dynamics<T, AD_Descriptor> *conductiveDynamics_,
        Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t) :
    conductiveDynamics(conductiveDynamics_), adiabaticDynamics(adiabaticDynamics_), iniVal(iniVal_)
{
    PLB_ASSERT(conductiveDynamics->getId() != adiabaticDynamics->getId());
    PLB_ASSERT(!util::isZero(Pr_t));
    C = FS_Descriptor<T>::cs2 * AD_Descriptor<T>::invCs2 / Pr_t;
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
AdvectionDiffusionFreeSurfaceCoupling3D<
    T, AD_Descriptor, FS_Descriptor>::~AdvectionDiffusionFreeSurfaceCoupling3D()
{
    delete conductiveDynamics;
    delete adiabaticDynamics;
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    AdvectionDiffusionFreeSurfaceCoupling3D(
        AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs) :
    conductiveDynamics(rhs.conductiveDynamics->clone()),
    adiabaticDynamics(rhs.adiabaticDynamics->clone()),
    iniVal(rhs.iniVal),
    C(rhs.C)
{ }

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>
    &AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::operator=(
        AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs)
{
    AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>(rhs).swap(*this);
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
void AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::swap(
    AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> &rhs)
{
    std::swap(conductiveDynamics, rhs.conductiveDynamics);
    std::swap(adiabaticDynamics, rhs.adiabaticDynamics);
    std::swap(iniVal, rhs.iniVal);
    std::swap(C, rhs.C);
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
void AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    using namespace freeSurfaceFlag;

    PLB_ASSERT(blocks.size() == 3);
    BlockLattice3D<T, AD_Descriptor> *tempLattice =
        dynamic_cast<BlockLattice3D<T, AD_Descriptor> *>(blocks[0]);
    BlockLattice3D<T, FS_Descriptor> *fluidLattice =
        dynamic_cast<BlockLattice3D<T, FS_Descriptor> *>(blocks[1]);
    ScalarField3D<int> *flags = dynamic_cast<ScalarField3D<int> *>(blocks[2]);
    PLB_ASSERT(tempLattice);
    PLB_ASSERT(fluidLattice);
    PLB_ASSERT(flags);

    plint conductiveId = conductiveDynamics->getId();
    plint adiabaticId = adiabaticDynamics->getId();
    BlockLattice3D<T, AD_Descriptor> tempBak(*tempLattice);

    Dot3D ofsFluid = computeRelativeDisplacement(*tempLattice, *fluidLattice);
    Dot3D ofsFlag = computeRelativeDisplacement(*tempLattice, *flags);

    const int velOffset = AD_Descriptor<T>::ExternalField::velocityBeginsAt;
    Array<T, 3> zero;
    zero.resetToZero();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint flag = flags->get(iX + ofsFlag.x, iY + ofsFlag.y, iZ + ofsFlag.z);
                bool tempIsAdiabatic = tempBak.get(iX, iY, iZ).getDynamics().getId() == adiabaticId;

                if (flag == fluid && tempIsAdiabatic) {
                    tempLattice->attributeDynamics(iX, iY, iZ, conductiveDynamics->clone());
                    Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY, iZ);
                    for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                        cell[iPop] = T();
                    }
                    plint numNeighbors = 0;
                    plint d = 1;
                    for (plint dx = -d; dx <= d; dx++) {
                        plint i = iX + dx;
                        for (plint dy = -d; dy <= d; dy++) {
                            plint j = iY + dy;
                            for (plint dz = -d; dz <= d; dz++) {
                                plint k = iZ + dz;
                                Cell<T, AD_Descriptor> const &nextCell = tempBak.get(i, j, k);
                                // plint nextFlag = flags->get(i+ofsFlag.x,j+ofsFlag.y,k+ofsFlag.z);
                                // if(!(dx==0 && dy==0 && dz==0) &&
                                // (nextCell.getDynamics().getId()==conductiveId) &&
                                //         (nextFlag!=wall)) {
                                if (!(dx == 0 && dy == 0 && dz == 0)
                                    && (nextCell.getDynamics().getId() == conductiveId)) {
                                    for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                                        cell[iPop] += nextCell[iPop];
                                    }
                                    ++numNeighbors;
                                }
                            }
                        }
                    }
                    for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                        if (numNeighbors == 0) {
                            cell[iPop] = AD_Descriptor<T>::t[iPop];
                        } else {
                            cell[iPop] /= numNeighbors;
                        }
                    }
                }

                if (flag == interface) {
                    tempLattice->attributeDynamics(iX, iY, iZ, conductiveDynamics->clone());
                    Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY, iZ);
                    iniCellAtEquilibrium(cell, iniVal, zero);
                    tempLattice->attributeDynamics(iX, iY, iZ, adiabaticDynamics->clone());
                }

                // Set the advection velocity.

                Cell<T, AD_Descriptor> &tempCell = tempLattice->get(iX, iY, iZ);
                Cell<T, FS_Descriptor> &fluidCell =
                    fluidLattice->get(iX + ofsFluid.x, iY + ofsFluid.y, iZ + ofsFluid.z);

                T *u = tempCell.getExternal(velOffset);
                if (isFullWet(flag)) {
                    Array<T, 3> velocity;
                    fluidCell.computeVelocity(velocity);
                    velocity.to_cArray(u);
                } else {
                    zero.to_cArray(u);
                }

                // Set the relaxation parameter for the advection-diffusion.

                if (tempCell.getDynamics().getId() != adiabaticId) {
                    T fluidOmega = fluidCell.getDynamics().getDynamicParameter(
                        dynamicParams::dynamicOmega, fluidCell);
                    if (!util::isZero(fluidOmega)) {
                        T fluidOmega0 = fluidCell.getDynamics().getOmega();
                        T tempOmega0 = tempCell.getDynamics().getOmega();
                        T tempOmega =
                            (T)1
                            / ((T)1 / tempOmega0 + C * ((T)1 / fluidOmega - (T)1 / fluidOmega0));
                        tempCell.getDynamics().setOmega(tempOmega);
                    }
                }
            }
        }
    }
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>
    *AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::clone() const
{
    return new AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>(*this);
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
void AdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;  // Temperature
    modified[1] = modif::nothing;        // FreeSurface Fluid
    modified[2] = modif::nothing;        // FreeSurface Flags
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    MaskedAdvectionDiffusionFreeSurfaceCoupling3D(
        Dynamics<T, AD_Descriptor> *conductiveDynamics_,
        Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t, int maskValue_) :
    conductiveDynamics(conductiveDynamics_),
    adiabaticDynamics(adiabaticDynamics_),
    iniVal(iniVal_),
    maskValue(maskValue_)
{
    PLB_ASSERT(conductiveDynamics->getId() != adiabaticDynamics->getId());
    PLB_ASSERT(!util::isZero(Pr_t));
    C = FS_Descriptor<T>::cs2 * AD_Descriptor<T>::invCs2 / Pr_t;
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
MaskedAdvectionDiffusionFreeSurfaceCoupling3D<
    T, AD_Descriptor, FS_Descriptor>::~MaskedAdvectionDiffusionFreeSurfaceCoupling3D()
{
    delete conductiveDynamics;
    delete adiabaticDynamics;
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    MaskedAdvectionDiffusionFreeSurfaceCoupling3D(
        MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs) :
    conductiveDynamics(rhs.conductiveDynamics->clone()),
    adiabaticDynamics(rhs.adiabaticDynamics->clone()),
    iniVal(rhs.iniVal),
    C(rhs.C),
    maskValue(rhs.maskValue)
{ }

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>
    &MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::operator=(
        MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs)
{
    MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>(rhs).swap(*this);
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
void MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::swap(
    MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> &rhs)
{
    std::swap(conductiveDynamics, rhs.conductiveDynamics);
    std::swap(adiabaticDynamics, rhs.adiabaticDynamics);
    std::swap(iniVal, rhs.iniVal);
    std::swap(C, rhs.C);
    std::swap(maskValue, rhs.maskValue);
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
void MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    using namespace freeSurfaceFlag;

    PLB_ASSERT(blocks.size() == 4);
    BlockLattice3D<T, AD_Descriptor> *tempLattice =
        dynamic_cast<BlockLattice3D<T, AD_Descriptor> *>(blocks[0]);
    BlockLattice3D<T, FS_Descriptor> *fluidLattice =
        dynamic_cast<BlockLattice3D<T, FS_Descriptor> *>(blocks[1]);
    ScalarField3D<int> *flags = dynamic_cast<ScalarField3D<int> *>(blocks[2]);
    ScalarField3D<int> *mask = dynamic_cast<ScalarField3D<int> *>(blocks[3]);
    PLB_ASSERT(tempLattice);
    PLB_ASSERT(fluidLattice);
    PLB_ASSERT(flags);
    PLB_ASSERT(mask);

    plint conductiveId = conductiveDynamics->getId();
    plint adiabaticId = adiabaticDynamics->getId();
    BlockLattice3D<T, AD_Descriptor> tempBak(*tempLattice);

    Dot3D ofsFluid = computeRelativeDisplacement(*tempLattice, *fluidLattice);
    Dot3D ofsFlag = computeRelativeDisplacement(*tempLattice, *flags);
    Dot3D ofsMask = computeRelativeDisplacement(*tempLattice, *mask);

    const int velOffset = AD_Descriptor<T>::ExternalField::velocityBeginsAt;
    Array<T, 3> zero;
    zero.resetToZero();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask->get(iX + ofsMask.x, iY + ofsMask.y, iZ + ofsMask.z) != maskValue) {
                    continue;
                }

                plint flag = flags->get(iX + ofsFlag.x, iY + ofsFlag.y, iZ + ofsFlag.z);
                bool tempIsAdiabatic = tempBak.get(iX, iY, iZ).getDynamics().getId() == adiabaticId;

                if (flag == fluid && tempIsAdiabatic) {
                    tempLattice->attributeDynamics(iX, iY, iZ, conductiveDynamics->clone());
                    Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY, iZ);
                    for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                        cell[iPop] = T();
                    }
                    plint numNeighbors = 0;
                    plint d = 1;
                    for (plint dx = -d; dx <= d; dx++) {
                        plint i = iX + dx;
                        for (plint dy = -d; dy <= d; dy++) {
                            plint j = iY + dy;
                            for (plint dz = -d; dz <= d; dz++) {
                                plint k = iZ + dz;
                                Cell<T, AD_Descriptor> const &nextCell = tempBak.get(i, j, k);
                                // plint nextFlag = flags->get(i+ofsFlag.x,j+ofsFlag.y,k+ofsFlag.z);
                                // if(!(dx==0 && dy==0 && dz==0) &&
                                // (nextCell.getDynamics().getId()==conductiveId) &&
                                //         (nextFlag!=wall)) {
                                //  TODO: Should the mask values be checked for the neighbors?
                                if (!(dx == 0 && dy == 0 && dz == 0)
                                    && (nextCell.getDynamics().getId() == conductiveId)) {
                                    for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                                        cell[iPop] += nextCell[iPop];
                                    }
                                    ++numNeighbors;
                                }
                            }
                        }
                    }
                    for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                        if (numNeighbors == 0) {
                            cell[iPop] = AD_Descriptor<T>::t[iPop];
                        } else {
                            cell[iPop] /= numNeighbors;
                        }
                    }
                }

                if (flag == interface) {
                    tempLattice->attributeDynamics(iX, iY, iZ, conductiveDynamics->clone());
                    Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY, iZ);
                    iniCellAtEquilibrium(cell, iniVal, zero);
                    tempLattice->attributeDynamics(iX, iY, iZ, adiabaticDynamics->clone());
                }

                // Set the advection velocity.

                Cell<T, AD_Descriptor> &tempCell = tempLattice->get(iX, iY, iZ);
                Cell<T, FS_Descriptor> &fluidCell =
                    fluidLattice->get(iX + ofsFluid.x, iY + ofsFluid.y, iZ + ofsFluid.z);

                T *u = tempCell.getExternal(velOffset);
                if (isFullWet(flag)) {
                    Array<T, 3> velocity;
                    fluidCell.computeVelocity(velocity);
                    velocity.to_cArray(u);
                } else {
                    zero.to_cArray(u);
                }

                // Set the relaxation parameter for the advection-diffusion.

                if (tempCell.getDynamics().getId() != adiabaticId) {
                    T fluidOmega = fluidCell.getDynamics().getDynamicParameter(
                        dynamicParams::dynamicOmega, fluidCell);
                    if (!util::isZero(fluidOmega)) {
                        T fluidOmega0 = fluidCell.getDynamics().getOmega();
                        T tempOmega0 = tempCell.getDynamics().getOmega();
                        T tempOmega =
                            (T)1
                            / ((T)1 / tempOmega0 + C * ((T)1 / fluidOmega - (T)1 / fluidOmega0));
                        tempCell.getDynamics().setOmega(tempOmega);
                    }
                }
            }
        }
    }
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>
    *MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::clone() const
{
    return new MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>(
        *this);
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
void MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;  // Temperature
    modified[1] = modif::nothing;        // FreeSurface Fluid
    modified[2] = modif::nothing;        // FreeSurface Flags
    modified[3] = modif::nothing;        // Mask
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D(
        Dynamics<T, AD_Descriptor> *conductiveDynamics_,
        Dynamics<T, AD_Descriptor> *adiabaticDynamics_, T iniVal_, T Pr_t, int maskValue_) :
    conductiveDynamics(conductiveDynamics_),
    adiabaticDynamics(adiabaticDynamics_),
    iniVal(iniVal_),
    maskValue(maskValue_)
{
    PLB_ASSERT(conductiveDynamics->getId() != adiabaticDynamics->getId());
    PLB_ASSERT(!util::isZero(Pr_t));
    C = FS_Descriptor<T>::cs2 * AD_Descriptor<T>::invCs2 / Pr_t;
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<
    T, AD_Descriptor, FS_Descriptor>::~N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D()
{
    delete conductiveDynamics;
    delete adiabaticDynamics;
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D(
        N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const
            &rhs) :
    conductiveDynamics(rhs.conductiveDynamics->clone()),
    adiabaticDynamics(rhs.adiabaticDynamics->clone()),
    iniVal(rhs.iniVal),
    C(rhs.C),
    maskValue(rhs.maskValue)
{ }

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>
    &N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::operator=(
        N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> const &rhs)
{
    N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>(rhs).swap(
        *this);
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
void N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::swap(
    N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor> &rhs)
{
    std::swap(conductiveDynamics, rhs.conductiveDynamics);
    std::swap(adiabaticDynamics, rhs.adiabaticDynamics);
    std::swap(iniVal, rhs.iniVal);
    std::swap(C, rhs.C);
    std::swap(maskValue, rhs.maskValue);
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
void N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    using namespace freeSurfaceFlag;

    PLB_ASSERT(blocks.size() == 4);
    BlockLattice3D<T, AD_Descriptor> *tempLattice =
        dynamic_cast<BlockLattice3D<T, AD_Descriptor> *>(blocks[0]);
    BlockLattice3D<T, FS_Descriptor> *fluidLattice =
        dynamic_cast<BlockLattice3D<T, FS_Descriptor> *>(blocks[1]);
    ScalarField3D<int> *flags = dynamic_cast<ScalarField3D<int> *>(blocks[2]);
    NTensorField3D<int> *mask = dynamic_cast<NTensorField3D<int> *>(blocks[3]);
    PLB_ASSERT(tempLattice);
    PLB_ASSERT(fluidLattice);
    PLB_ASSERT(flags);
    PLB_ASSERT(mask);

    plint conductiveId = conductiveDynamics->getId();
    plint adiabaticId = adiabaticDynamics->getId();
    BlockLattice3D<T, AD_Descriptor> tempBak(*tempLattice);

    Dot3D ofsFluid = computeRelativeDisplacement(*tempLattice, *fluidLattice);
    Dot3D ofsFlag = computeRelativeDisplacement(*tempLattice, *flags);
    Dot3D ofsMask = computeRelativeDisplacement(*tempLattice, *mask);

    const int velOffset = AD_Descriptor<T>::ExternalField::velocityBeginsAt;
    Array<T, 3> zero;
    zero.resetToZero();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (*mask->get(iX + ofsMask.x, iY + ofsMask.y, iZ + ofsMask.z) != maskValue) {
                    continue;
                }

                plint flag = flags->get(iX + ofsFlag.x, iY + ofsFlag.y, iZ + ofsFlag.z);
                bool tempIsAdiabatic = tempBak.get(iX, iY, iZ).getDynamics().getId() == adiabaticId;

                if (flag == fluid && tempIsAdiabatic) {
                    tempLattice->attributeDynamics(iX, iY, iZ, conductiveDynamics->clone());
                    Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY, iZ);
                    for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                        cell[iPop] = T();
                    }
                    plint numNeighbors = 0;
                    plint d = 1;
                    for (plint dx = -d; dx <= d; dx++) {
                        plint i = iX + dx;
                        for (plint dy = -d; dy <= d; dy++) {
                            plint j = iY + dy;
                            for (plint dz = -d; dz <= d; dz++) {
                                plint k = iZ + dz;
                                Cell<T, AD_Descriptor> const &nextCell = tempBak.get(i, j, k);
                                // plint nextFlag = flags->get(i+ofsFlag.x,j+ofsFlag.y,k+ofsFlag.z);
                                // if(!(dx==0 && dy==0 && dz==0) &&
                                // (nextCell.getDynamics().getId()==conductiveId) &&
                                //         (nextFlag!=wall)) {
                                //  TODO: Should the mask values be checked for the neighbors?
                                if (!(dx == 0 && dy == 0 && dz == 0)
                                    && (nextCell.getDynamics().getId() == conductiveId)) {
                                    for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                                        cell[iPop] += nextCell[iPop];
                                    }
                                    ++numNeighbors;
                                }
                            }
                        }
                    }
                    for (plint iPop = 0; iPop < AD_Descriptor<T>::q; ++iPop) {
                        if (numNeighbors == 0) {
                            cell[iPop] = AD_Descriptor<T>::t[iPop];
                        } else {
                            cell[iPop] /= numNeighbors;
                        }
                    }
                }

                if (flag == interface) {
                    tempLattice->attributeDynamics(iX, iY, iZ, conductiveDynamics->clone());
                    Cell<T, AD_Descriptor> &cell = tempLattice->get(iX, iY, iZ);
                    iniCellAtEquilibrium(cell, iniVal, zero);
                    tempLattice->attributeDynamics(iX, iY, iZ, adiabaticDynamics->clone());
                }

                // Set the advection velocity.

                Cell<T, AD_Descriptor> &tempCell = tempLattice->get(iX, iY, iZ);
                Cell<T, FS_Descriptor> &fluidCell =
                    fluidLattice->get(iX + ofsFluid.x, iY + ofsFluid.y, iZ + ofsFluid.z);

                T *u = tempCell.getExternal(velOffset);
                if (isFullWet(flag)) {
                    Array<T, 3> velocity;
                    fluidCell.computeVelocity(velocity);
                    velocity.to_cArray(u);
                } else {
                    zero.to_cArray(u);
                }

                // Set the relaxation parameter for the advection-diffusion.

                if (tempCell.getDynamics().getId() != adiabaticId) {
                    T fluidOmega = fluidCell.getDynamics().getDynamicParameter(
                        dynamicParams::dynamicOmega, fluidCell);
                    if (!util::isZero(fluidOmega)) {
                        T fluidOmega0 = fluidCell.getDynamics().getOmega();
                        T tempOmega0 = tempCell.getDynamics().getOmega();
                        T tempOmega =
                            (T)1
                            / ((T)1 / tempOmega0 + C * ((T)1 / fluidOmega - (T)1 / fluidOmega0));
                        tempCell.getDynamics().setOmega(tempOmega);
                    }
                }
            }
        }
    }
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>
    *N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::clone() const
{
    return new N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>(
        *this);
}

template <
    typename T, template <typename U1> class AD_Descriptor,
    template <typename U2> class FS_Descriptor>
void N_MaskedAdvectionDiffusionFreeSurfaceCoupling3D<T, AD_Descriptor, FS_Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::dataStructure;  // Temperature
    modified[1] = modif::nothing;        // FreeSurface Fluid
    modified[2] = modif::nothing;        // FreeSurface Flags
    modified[3] = modif::nothing;        // Mask
}

template <typename T>
AdvectionDiffusionFd3D<T>::AdvectionDiffusionFd3D(T d_, bool upwind_) : d(d_), upwind(upwind_)
{ }

template <typename T>
void AdvectionDiffusionFd3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 5);
    ScalarField3D<T> *phi_t = dynamic_cast<ScalarField3D<T> *>(fields[0]);
    ScalarField3D<T> *phi_tp1 = dynamic_cast<ScalarField3D<T> *>(fields[1]);
    ScalarField3D<T> *result = dynamic_cast<ScalarField3D<T> *>(fields[2]);
    TensorField3D<T, 3> *uField = dynamic_cast<TensorField3D<T, 3> *>(fields[3]);
    ScalarField3D<T> *Q = dynamic_cast<ScalarField3D<T> *>(fields[4]);
    Dot3D ofs1 = computeRelativeDisplacement(*phi_t, *phi_tp1);
    Dot3D ofs2 = computeRelativeDisplacement(*phi_t, *result);
    Dot3D ofs3 = computeRelativeDisplacement(*phi_t, *uField);
    Dot3D ofs4 = computeRelativeDisplacement(*phi_t, *Q);

    if (upwind) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    T phiC = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiE = phi_tp1->get(iX + 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiW = phi_tp1->get(iX - 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiN = phi_tp1->get(iX + ofs1.x, iY + 1 + ofs1.y, iZ + ofs1.z);
                    T phiS = phi_tp1->get(iX + ofs1.x, iY - 1 + ofs1.y, iZ + ofs1.z);
                    T phiT = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + 1 + ofs1.z);
                    T phiB = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ - 1 + ofs1.z);

                    Array<T, 3> const &u = uField->get(iX + ofs3.x, iY + ofs3.y, iZ + ofs3.z);
                    T advection =
                        u[0]
                            * (util::greaterThan(u[0], (T)0)
                                   ? (phiC - phiW)
                                   : (util::lessThan(u[0], (T)0) ? (phiE - phiC)
                                                                 : (T)0.5 * (phiE - phiW)))
                        + u[1]
                              * (util::greaterThan(u[1], (T)0)
                                     ? (phiC - phiS)
                                     : (util::lessThan(u[1], (T)0) ? (phiN - phiC)
                                                                   : (T)0.5 * (phiN - phiS)))
                        + u[2]
                              * (util::greaterThan(u[2], (T)0)
                                     ? (phiC - phiB)
                                     : (util::lessThan(u[2], (T)0) ? (phiT - phiC)
                                                                   : (T)0.5 * (phiT - phiB)));

                    T diffusion = d * (phiE + phiW + phiN + phiS + phiT + phiB - (T)6 * phiC);

                    result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) =
                        phi_t->get(iX, iY, iZ) + diffusion - advection
                        + Q->get(iX + ofs4.x, iY + ofs4.y, iZ + ofs4.z);
                }
            }
        }
    } else {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    T phiC = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiE = phi_tp1->get(iX + 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiW = phi_tp1->get(iX - 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiN = phi_tp1->get(iX + ofs1.x, iY + 1 + ofs1.y, iZ + ofs1.z);
                    T phiS = phi_tp1->get(iX + ofs1.x, iY - 1 + ofs1.y, iZ + ofs1.z);
                    T phiT = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + 1 + ofs1.z);
                    T phiB = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ - 1 + ofs1.z);

                    Array<T, 3> const &u = uField->get(iX + ofs3.x, iY + ofs3.y, iZ + ofs3.z);
                    T advection =
                        (T)0.5
                        * (u[0] * (phiE - phiW) + u[1] * (phiN - phiS) + u[2] * (phiT - phiB));

                    T diffusion = d * (phiE + phiW + phiN + phiS + phiT + phiB - (T)6 * phiC);

                    result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) =
                        phi_t->get(iX, iY, iZ) + diffusion - advection
                        + Q->get(iX + ofs4.x, iY + ofs4.y, iZ + ofs4.z);
                }
            }
        }
    }
}

template <typename T>
AdvectionDiffusionFd3D<T> *AdvectionDiffusionFd3D<T>::clone() const
{
    return new AdvectionDiffusionFd3D<T>(*this);
}

template <typename T>
void AdvectionDiffusionFd3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // phi_t
    modified[1] = modif::nothing;          // phi_tp1
    modified[2] = modif::staticVariables;  // result
    modified[3] = modif::nothing;          // u
    modified[4] = modif::nothing;          // Q
}

template <typename T>
VariableDiffusivityAdvectionDiffusionFd3D<T>::VariableDiffusivityAdvectionDiffusionFd3D(
    bool upwind_) :
    upwind(upwind_)
{ }

template <typename T>
void VariableDiffusivityAdvectionDiffusionFd3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 6);
    ScalarField3D<T> *phi_t = dynamic_cast<ScalarField3D<T> *>(fields[0]);
    ScalarField3D<T> *phi_tp1 = dynamic_cast<ScalarField3D<T> *>(fields[1]);
    ScalarField3D<T> *result = dynamic_cast<ScalarField3D<T> *>(fields[2]);
    TensorField3D<T, 3> *uField = dynamic_cast<TensorField3D<T, 3> *>(fields[3]);
    ScalarField3D<T> *Q = dynamic_cast<ScalarField3D<T> *>(fields[4]);
    ScalarField3D<T> *D = dynamic_cast<ScalarField3D<T> *>(fields[5]);
    Dot3D ofs1 = computeRelativeDisplacement(*phi_t, *phi_tp1);
    Dot3D ofs2 = computeRelativeDisplacement(*phi_t, *result);
    Dot3D ofs3 = computeRelativeDisplacement(*phi_t, *uField);
    Dot3D ofs4 = computeRelativeDisplacement(*phi_t, *Q);
    Dot3D ofs5 = computeRelativeDisplacement(*phi_t, *D);

    if (upwind) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    T phiC = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiE = phi_tp1->get(iX + 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiW = phi_tp1->get(iX - 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiN = phi_tp1->get(iX + ofs1.x, iY + 1 + ofs1.y, iZ + ofs1.z);
                    T phiS = phi_tp1->get(iX + ofs1.x, iY - 1 + ofs1.y, iZ + ofs1.z);
                    T phiT = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + 1 + ofs1.z);
                    T phiB = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ - 1 + ofs1.z);

                    Array<T, 3> const &u = uField->get(iX + ofs3.x, iY + ofs3.y, iZ + ofs3.z);
                    T advection =
                        u[0]
                            * (util::greaterThan(u[0], (T)0)
                                   ? (phiC - phiW)
                                   : (util::lessThan(u[0], (T)0) ? (phiE - phiC)
                                                                 : (T)0.5 * (phiE - phiW)))
                        + u[1]
                              * (util::greaterThan(u[1], (T)0)
                                     ? (phiC - phiS)
                                     : (util::lessThan(u[1], (T)0) ? (phiN - phiC)
                                                                   : (T)0.5 * (phiN - phiS)))
                        + u[2]
                              * (util::greaterThan(u[2], (T)0)
                                     ? (phiC - phiB)
                                     : (util::lessThan(u[2], (T)0) ? (phiT - phiC)
                                                                   : (T)0.5 * (phiT - phiB)));

                    T DC = D->get(iX + ofs5.x, iY + ofs5.y, iZ + ofs5.z);
                    T DE = (T)0.5 * (DC + D->get(iX + 1 + ofs5.x, iY + ofs5.y, iZ + ofs5.z));
                    T DW = (T)0.5 * (DC + D->get(iX - 1 + ofs5.x, iY + ofs5.y, iZ + ofs5.z));
                    T DN = (T)0.5 * (DC + D->get(iX + ofs5.x, iY + 1 + ofs5.y, iZ + ofs5.z));
                    T DS = (T)0.5 * (DC + D->get(iX + ofs5.x, iY - 1 + ofs5.y, iZ + ofs5.z));
                    T DT = (T)0.5 * (DC + D->get(iX + ofs5.x, iY + ofs5.y, iZ + 1 + ofs5.z));
                    T DB = (T)0.5 * (DC + D->get(iX + ofs5.x, iY + ofs5.y, iZ - 1 + ofs5.z));

                    // The diffusion term is implemented as: div(D * grad(T)).
                    T diffusion = DE * (phiE - phiC) - DW * (phiC - phiW) + DN * (phiN - phiC)
                                  - DS * (phiC - phiS) + DT * (phiT - phiC) - DB * (phiC - phiB);

                    result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) =
                        phi_t->get(iX, iY, iZ) + diffusion - advection
                        + Q->get(iX + ofs4.x, iY + ofs4.y, iZ + ofs4.z);
                }
            }
        }
    } else {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    T phiC = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiE = phi_tp1->get(iX + 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiW = phi_tp1->get(iX - 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    T phiN = phi_tp1->get(iX + ofs1.x, iY + 1 + ofs1.y, iZ + ofs1.z);
                    T phiS = phi_tp1->get(iX + ofs1.x, iY - 1 + ofs1.y, iZ + ofs1.z);
                    T phiT = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + 1 + ofs1.z);
                    T phiB = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ - 1 + ofs1.z);

                    Array<T, 3> const &u = uField->get(iX + ofs3.x, iY + ofs3.y, iZ + ofs3.z);
                    T advection =
                        (T)0.5
                        * (u[0] * (phiE - phiW) + u[1] * (phiN - phiS) + u[2] * (phiT - phiB));

                    T DC = D->get(iX + ofs5.x, iY + ofs5.y, iZ + ofs5.z);
                    T DE = (T)0.5 * (DC + D->get(iX + 1 + ofs5.x, iY + ofs5.y, iZ + ofs5.z));
                    T DW = (T)0.5 * (DC + D->get(iX - 1 + ofs5.x, iY + ofs5.y, iZ + ofs5.z));
                    T DN = (T)0.5 * (DC + D->get(iX + ofs5.x, iY + 1 + ofs5.y, iZ + ofs5.z));
                    T DS = (T)0.5 * (DC + D->get(iX + ofs5.x, iY - 1 + ofs5.y, iZ + ofs5.z));
                    T DT = (T)0.5 * (DC + D->get(iX + ofs5.x, iY + ofs5.y, iZ + 1 + ofs5.z));
                    T DB = (T)0.5 * (DC + D->get(iX + ofs5.x, iY + ofs5.y, iZ - 1 + ofs5.z));

                    // The diffusion term is implemented as: div(D * grad(T)).
                    T diffusion = DE * (phiE - phiC) - DW * (phiC - phiW) + DN * (phiN - phiC)
                                  - DS * (phiC - phiS) + DT * (phiT - phiC) - DB * (phiC - phiB);

                    result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) =
                        phi_t->get(iX, iY, iZ) + diffusion - advection
                        + Q->get(iX + ofs4.x, iY + ofs4.y, iZ + ofs4.z);
                }
            }
        }
    }
}

template <typename T>
VariableDiffusivityAdvectionDiffusionFd3D<T> *VariableDiffusivityAdvectionDiffusionFd3D<T>::clone()
    const
{
    return new VariableDiffusivityAdvectionDiffusionFd3D<T>(*this);
}

template <typename T>
void VariableDiffusivityAdvectionDiffusionFd3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // phi_t
    modified[1] = modif::nothing;          // phi_tp1
    modified[2] = modif::staticVariables;  // result
    modified[3] = modif::nothing;          // u
    modified[4] = modif::nothing;          // Q
    modified[5] = modif::nothing;          // D
}

template <typename T, class ADCellT>
MaskedVariableDiffusivityAdvectionDiffusionFd3D<T, ADCellT>::
    MaskedVariableDiffusivityAdvectionDiffusionFd3D(
        ADCellT adCellT_, T defaultValue_, Box3D const &fullClosedDomain_, bool xPeriodic_,
        bool yPeriodic_, bool zPeriodic_, bool upwind_) :
    adCellT(adCellT_),
    defaultValue(defaultValue_),
    fullExtendedDomain(fullClosedDomain_),
    upwind(upwind_)
{
    if (xPeriodic_) {
        fullExtendedDomain.x0--;
        fullExtendedDomain.x1++;
    }
    if (yPeriodic_) {
        fullExtendedDomain.y0--;
        fullExtendedDomain.y1++;
    }
    if (zPeriodic_) {
        fullExtendedDomain.z0--;
        fullExtendedDomain.z1++;
    }
}

template <typename T, class ADCellT>
void MaskedVariableDiffusivityAdvectionDiffusionFd3D<T, ADCellT>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 7);
    ScalarField3D<T> *phi_t = dynamic_cast<ScalarField3D<T> *>(fields[0]);
    ScalarField3D<T> *phi_tp1 = dynamic_cast<ScalarField3D<T> *>(fields[1]);
    ScalarField3D<T> *result = dynamic_cast<ScalarField3D<T> *>(fields[2]);
    TensorField3D<T, 3> *uField = dynamic_cast<TensorField3D<T, 3> *>(fields[3]);
    ScalarField3D<T> *Q = dynamic_cast<ScalarField3D<T> *>(fields[4]);
    ScalarField3D<T> *D = dynamic_cast<ScalarField3D<T> *>(fields[5]);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int> *>(fields[6]);
    Dot3D ofs1 = computeRelativeDisplacement(*phi_t, *phi_tp1);
    Dot3D ofs2 = computeRelativeDisplacement(*phi_t, *result);
    Dot3D ofs3 = computeRelativeDisplacement(*phi_t, *uField);
    Dot3D ofs4 = computeRelativeDisplacement(*phi_t, *Q);
    Dot3D ofs5 = computeRelativeDisplacement(*phi_t, *D);
    Dot3D ofs6 = computeRelativeDisplacement(*phi_t, *flag);

    Dot3D absOfs = phi_t->getLocation();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                int flagC = flag->get(iX + ofs6.x, iY + ofs6.y, iZ + ofs6.z);
                T phi = (T)0;
                T dphi = (T)0;
                if (adCellT.isSingular(flagC)) {
                    result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) = defaultValue;
                } else if (adCellT.isFixed(flagC, phi)) {
                    result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) = phi;
                } else if (adCellT.isZeroFlux(flagC)) {
                    T avePhi = (T)0;
                    plint n = 0;
                    for (plint dx = -1; dx <= 1; ++dx) {
                        plint nextX = iX + dx;
                        for (plint dy = -1; dy <= 1; ++dy) {
                            plint nextY = iY + dy;
                            for (plint dz = -1; dz <= 1; ++dz) {
                                plint nextZ = iZ + dz;
                                if (contained(
                                        nextX + absOfs.x, nextY + absOfs.y, nextZ + absOfs.z,
                                        fullExtendedDomain)) {
                                    int nextFlag =
                                        flag->get(nextX + ofs6.x, nextY + ofs6.y, nextZ + ofs6.z);
                                    if (!adCellT.isSingular(nextFlag)
                                        && !adCellT.isZeroFlux(nextFlag)
                                        && !adCellT.isFixedFlux(nextFlag, phi, dphi))
                                    {
                                        avePhi += phi_tp1->get(
                                            nextX + ofs1.x, nextY + ofs1.y, nextZ + ofs1.z);
                                        n++;
                                    }
                                }
                            }
                        }
                    }
                    if (n != 0) {
                        avePhi /= (T)n;
                        result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) = avePhi;
                    } else {
                        result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) = defaultValue;
                    }
                } else if (adCellT.isFixedFlux(flagC, phi, dphi)) {
                    T avePhi = (T)0;
                    plint n = 0;
                    for (plint dx = -1; dx <= 1; ++dx) {
                        plint nextX = iX + dx;
                        for (plint dy = -1; dy <= 1; ++dy) {
                            plint nextY = iY + dy;
                            for (plint dz = -1; dz <= 1; ++dz) {
                                plint nextZ = iZ + dz;
                                if (contained(
                                        nextX + absOfs.x, nextY + absOfs.y, nextZ + absOfs.z,
                                        fullExtendedDomain)) {
                                    int nextFlag =
                                        flag->get(nextX + ofs6.x, nextY + ofs6.y, nextZ + ofs6.z);
                                    if (!adCellT.isSingular(nextFlag)
                                        && !adCellT.isZeroFlux(nextFlag)
                                        && !adCellT.isFixedFlux(nextFlag, phi, dphi))
                                    {
                                        phi = phi_tp1->get(
                                            nextX + ofs1.x, nextY + ofs1.y, nextZ + ofs1.z);
                                        dphi = (T)0;
                                        (void)adCellT.isFixedFlux(flagC, phi, dphi);
                                        avePhi += phi + dphi;
                                        n++;
                                    }
                                }
                            }
                        }
                    }
                    if (n != 0) {
                        avePhi /= (T)n;
                        result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) = avePhi;
                    } else {
                        result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) = defaultValue;
                    }
                } else {
                    int flagE = flag->get(iX + 1 + ofs6.x, iY + ofs6.y, iZ + ofs6.z);
                    int flagW = flag->get(iX - 1 + ofs6.x, iY + ofs6.y, iZ + ofs6.z);
                    int flagN = flag->get(iX + ofs6.x, iY + 1 + ofs6.y, iZ + ofs6.z);
                    int flagS = flag->get(iX + ofs6.x, iY - 1 + ofs6.y, iZ + ofs6.z);
                    int flagT = flag->get(iX + ofs6.x, iY + ofs6.y, iZ + 1 + ofs6.z);
                    int flagB = flag->get(iX + ofs6.x, iY + ofs6.y, iZ - 1 + ofs6.z);

                    T phiC = phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + ofs1.z);
                    phi = dphi = (T)0;

                    T phiE = adCellT.isFixed(flagE, phi)
                                 ? phi
                                 : (adCellT.isSingular(flagE) || adCellT.isZeroFlux(flagE)
                                            || adCellT.isFixedFlux(flagE, phiC, dphi)
                                        ? phiC + dphi
                                        : phi_tp1->get(iX + 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z));
                    phi = dphi = (T)0;
                    T phiW = adCellT.isFixed(flagW, phi)
                                 ? phi
                                 : (adCellT.isSingular(flagW) || adCellT.isZeroFlux(flagW)
                                            || adCellT.isFixedFlux(flagW, phiC, dphi)
                                        ? phiC + dphi
                                        : phi_tp1->get(iX - 1 + ofs1.x, iY + ofs1.y, iZ + ofs1.z));
                    phi = dphi = (T)0;
                    T phiN = adCellT.isFixed(flagN, phi)
                                 ? phi
                                 : (adCellT.isSingular(flagN) || adCellT.isZeroFlux(flagN)
                                            || adCellT.isFixedFlux(flagN, phiC, dphi)
                                        ? phiC + dphi
                                        : phi_tp1->get(iX + ofs1.x, iY + 1 + ofs1.y, iZ + ofs1.z));
                    phi = dphi = (T)0;
                    T phiS = adCellT.isFixed(flagS, phi)
                                 ? phi
                                 : (adCellT.isSingular(flagS) || adCellT.isZeroFlux(flagS)
                                            || adCellT.isFixedFlux(flagS, phiC, dphi)
                                        ? phiC + dphi
                                        : phi_tp1->get(iX + ofs1.x, iY - 1 + ofs1.y, iZ + ofs1.z));
                    phi = dphi = (T)0;
                    T phiT = adCellT.isFixed(flagT, phi)
                                 ? phi
                                 : (adCellT.isSingular(flagT) || adCellT.isZeroFlux(flagT)
                                            || adCellT.isFixedFlux(flagT, phiC, dphi)
                                        ? phiC + dphi
                                        : phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ + 1 + ofs1.z));
                    phi = dphi = (T)0;
                    T phiB = adCellT.isFixed(flagB, phi)
                                 ? phi
                                 : (adCellT.isSingular(flagB) || adCellT.isZeroFlux(flagB)
                                            || adCellT.isFixedFlux(flagB, phiC, dphi)
                                        ? phiC + dphi
                                        : phi_tp1->get(iX + ofs1.x, iY + ofs1.y, iZ - 1 + ofs1.z));
                    phi = dphi = (T)0;

                    Array<T, 3> const &u = uField->get(iX + ofs3.x, iY + ofs3.y, iZ + ofs3.z);
                    T advection = (T)0;
                    if (upwind) {
                        advection =
                            u[0]
                                * (util::greaterThan(u[0], (T)0)
                                       ? (phiC - phiW)
                                       : (util::lessThan(u[0], (T)0) ? (phiE - phiC)
                                                                     : (T)0.5 * (phiE - phiW)))
                            + u[1]
                                  * (util::greaterThan(u[1], (T)0)
                                         ? (phiC - phiS)
                                         : (util::lessThan(u[1], (T)0) ? (phiN - phiC)
                                                                       : (T)0.5 * (phiN - phiS)))
                            + u[2]
                                  * (util::greaterThan(u[2], (T)0)
                                         ? (phiC - phiB)
                                         : (util::lessThan(u[2], (T)0) ? (phiT - phiC)
                                                                       : (T)0.5 * (phiT - phiB)));
                    } else {
                        advection =
                            (T)0.5
                            * (u[0] * (phiE - phiW) + u[1] * (phiN - phiS) + u[2] * (phiT - phiB));
                    }

                    T DC = D->get(iX + ofs5.x, iY + ofs5.y, iZ + ofs5.z);
                    T DE = (T)0.5 * (DC + D->get(iX + 1 + ofs5.x, iY + ofs5.y, iZ + ofs5.z));
                    T DW = (T)0.5 * (DC + D->get(iX - 1 + ofs5.x, iY + ofs5.y, iZ + ofs5.z));
                    T DN = (T)0.5 * (DC + D->get(iX + ofs5.x, iY + 1 + ofs5.y, iZ + ofs5.z));
                    T DS = (T)0.5 * (DC + D->get(iX + ofs5.x, iY - 1 + ofs5.y, iZ + ofs5.z));
                    T DT = (T)0.5 * (DC + D->get(iX + ofs5.x, iY + ofs5.y, iZ + 1 + ofs5.z));
                    T DB = (T)0.5 * (DC + D->get(iX + ofs5.x, iY + ofs5.y, iZ - 1 + ofs5.z));

                    // The diffusion term is implemented as: div(D * grad(T)).
                    T diffusion = DE * (phiE - phiC) - DW * (phiC - phiW) + DN * (phiN - phiC)
                                  - DS * (phiC - phiS) + DT * (phiT - phiC) - DB * (phiC - phiB);

                    result->get(iX + ofs2.x, iY + ofs2.y, iZ + ofs2.z) =
                        phi_t->get(iX, iY, iZ) + diffusion - advection
                        + Q->get(iX + ofs4.x, iY + ofs4.y, iZ + ofs4.z);
                }
            }
        }
    }
}

template <typename T, class ADCellT>
MaskedVariableDiffusivityAdvectionDiffusionFd3D<T, ADCellT>
    *MaskedVariableDiffusivityAdvectionDiffusionFd3D<T, ADCellT>::clone() const
{
    return new MaskedVariableDiffusivityAdvectionDiffusionFd3D<T, ADCellT>(*this);
}

template <typename T, class ADCellT>
void MaskedVariableDiffusivityAdvectionDiffusionFd3D<T, ADCellT>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // phi_t
    modified[1] = modif::nothing;          // phi_tp1
    modified[2] = modif::staticVariables;  // result
    modified[3] = modif::nothing;          // u
    modified[4] = modif::nothing;          // Q
    modified[5] = modif::nothing;          // D
    modified[6] = modif::nothing;          // flag
}

template <typename T, template <typename U> class FluidDescriptor>
TurbulentDiffusivityFunctional3D<T, FluidDescriptor>::TurbulentDiffusivityFunctional3D(
    T D0_, T Pr_t) :
    D0(D0_)
{
    PLB_ASSERT(!util::isZero(Pr_t));
    C = FluidDescriptor<T>::cs2 / Pr_t;
}

template <typename T, template <typename U> class FluidDescriptor>
void TurbulentDiffusivityFunctional3D<T, FluidDescriptor>::process(
    Box3D domain, BlockLattice3D<T, FluidDescriptor> &fluid, ScalarField3D<T> &diffusivity)
{
    Dot3D offset = computeRelativeDisplacement(fluid, diffusivity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T D = D0;
                Cell<T, FluidDescriptor> &cell = fluid.get(iX, iY, iZ);
                T omega = cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                if (!util::isZero(omega)) {
                    T omega0 = cell.getDynamics().getOmega();
                    D += C * ((T)1 / omega - (T)1 / omega0);
                }
                diffusivity.get(iX + offset.x, iY + offset.y, iZ + offset.z) = D;
            }
        }
    }
}

template <typename T, template <typename U> class FluidDescriptor>
TurbulentDiffusivityFunctional3D<T, FluidDescriptor>
    *TurbulentDiffusivityFunctional3D<T, FluidDescriptor>::clone() const
{
    return new TurbulentDiffusivityFunctional3D<T, FluidDescriptor>(*this);
}

template <typename T, template <typename U> class FluidDescriptor>
void TurbulentDiffusivityFunctional3D<T, FluidDescriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::staticVariables;  // Diffusivity.
}

template <typename T, template <typename U> class FluidDescriptor>
void computeTurbulentDiffusivity(
    T D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiScalarField3D<T> &diffusivity, Box3D domain)
{
    applyProcessingFunctional(
        new TurbulentDiffusivityFunctional3D<T, FluidDescriptor>(D0, Pr_t), domain, fluid,
        diffusivity);
}

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivity(
    T D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > diffusivity =
        generateMultiScalarField<T>(fluid, domain);
    computeTurbulentDiffusivity(D0, Pr_t, fluid, *diffusivity, domain);
    return diffusivity;
}

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivity(
    T D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid)
{
    return computeTurbulentDiffusivity(D0, Pr_t, fluid, fluid.getBoundingBox());
}

template <typename T, template <typename U> class FluidDescriptor>
TurbulentDiffusivityFromFlagFunctional3D<T, FluidDescriptor>::
    TurbulentDiffusivityFromFlagFunctional3D(std::vector<T> const &D0_, T Pr_t) :
    D0(D0_)
{
    PLB_ASSERT(!util::isZero(Pr_t));
    C = FluidDescriptor<T>::cs2 / Pr_t;
}

template <typename T, template <typename U> class FluidDescriptor>
void TurbulentDiffusivityFromFlagFunctional3D<T, FluidDescriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);

    BlockLattice3D<T, FluidDescriptor> *fluid =
        dynamic_cast<BlockLattice3D<T, FluidDescriptor> *>(blocks[0]);
    PLB_ASSERT(fluid);

    ScalarField3D<int> *flags = dynamic_cast<ScalarField3D<int> *>(blocks[1]);
    PLB_ASSERT(flags);

    ScalarField3D<T> *diffusivity = dynamic_cast<ScalarField3D<T> *>(blocks[2]);
    PLB_ASSERT(diffusivity);

    Dot3D ofsF = computeRelativeDisplacement(*fluid, *flags);
    Dot3D ofsD = computeRelativeDisplacement(*fluid, *diffusivity);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                int flag = flags->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z);
                if (flag >= 0) {
                    T D = D0[flag];
                    Cell<T, FluidDescriptor> &cell = fluid->get(iX, iY, iZ);
                    T omega =
                        cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                    if (!util::isZero(omega)) {
                        T omega0 = cell.getDynamics().getOmega();
                        D += C * ((T)1 / omega - (T)1 / omega0);
                    }
                    diffusivity->get(iX + ofsD.x, iY + ofsD.y, iZ + ofsD.z) = D;
                }
            }
        }
    }
}

template <typename T, template <typename U> class FluidDescriptor>
TurbulentDiffusivityFromFlagFunctional3D<T, FluidDescriptor>
    *TurbulentDiffusivityFromFlagFunctional3D<T, FluidDescriptor>::clone() const
{
    return new TurbulentDiffusivityFromFlagFunctional3D<T, FluidDescriptor>(*this);
}

template <typename T, template <typename U> class FluidDescriptor>
void TurbulentDiffusivityFromFlagFunctional3D<T, FluidDescriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::nothing;          // Flags.
    modified[2] = modif::staticVariables;  // Diffusivity.
}

template <typename T, template <typename U> class FluidDescriptor>
void computeTurbulentDiffusivityFromFlag(
    std::vector<T> const &D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiScalarField3D<int> &flags, MultiScalarField3D<T> &diffusivity, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&fluid);
    args.push_back(&flags);
    args.push_back(&diffusivity);
    applyProcessingFunctional(
        new TurbulentDiffusivityFromFlagFunctional3D<T, FluidDescriptor>(D0, Pr_t), domain, args);
}

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivityFromFlag(
    std::vector<T> const &D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiScalarField3D<int> &flags, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > diffusivity =
        generateMultiScalarField<T>(fluid, domain);
    computeTurbulentDiffusivityFromFlag(D0, Pr_t, fluid, flags, *diffusivity, domain);
    return diffusivity;
}

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivityFromFlag(
    std::vector<T> const &D0, T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid,
    MultiScalarField3D<int> &flags)
{
    return computeTurbulentDiffusivityFromFlag(D0, Pr_t, fluid, flags, fluid.getBoundingBox());
}

template <typename T, template <typename U> class FluidDescriptor>
TurbulentDiffusivityFromScalarFunctional3D<
    T, FluidDescriptor>::TurbulentDiffusivityFromScalarFunctional3D(T Pr_t)
{
    PLB_ASSERT(!util::isZero(Pr_t));
    C = FluidDescriptor<T>::cs2 / Pr_t;
}

template <typename T, template <typename U> class FluidDescriptor>
void TurbulentDiffusivityFromScalarFunctional3D<T, FluidDescriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);

    BlockLattice3D<T, FluidDescriptor> *fluid =
        dynamic_cast<BlockLattice3D<T, FluidDescriptor> *>(blocks[0]);
    PLB_ASSERT(fluid);

    ScalarField3D<T> *D0 = dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    PLB_ASSERT(D0);

    ScalarField3D<T> *diffusivity = dynamic_cast<ScalarField3D<T> *>(blocks[2]);
    PLB_ASSERT(diffusivity);

    Dot3D ofsD0 = computeRelativeDisplacement(*fluid, *D0);
    Dot3D ofsD = computeRelativeDisplacement(*fluid, *diffusivity);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T D = D0->get(iX + ofsD0.x, iY + ofsD0.y, iZ + ofsD0.z);
                Cell<T, FluidDescriptor> &cell = fluid->get(iX, iY, iZ);
                T omega = cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                if (!util::isZero(omega)) {
                    T omega0 = cell.getDynamics().getOmega();
                    D += C * ((T)1 / omega - (T)1 / omega0);
                }
                diffusivity->get(iX + ofsD.x, iY + ofsD.y, iZ + ofsD.z) = D;
            }
        }
    }
}

template <typename T, template <typename U> class FluidDescriptor>
TurbulentDiffusivityFromScalarFunctional3D<T, FluidDescriptor>
    *TurbulentDiffusivityFromScalarFunctional3D<T, FluidDescriptor>::clone() const
{
    return new TurbulentDiffusivityFromScalarFunctional3D<T, FluidDescriptor>(*this);
}

template <typename T, template <typename U> class FluidDescriptor>
void TurbulentDiffusivityFromScalarFunctional3D<T, FluidDescriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Fluid.
    modified[1] = modif::nothing;          // Base diffusivity.
    modified[2] = modif::staticVariables;  // Diffusivity.
}

template <typename T, template <typename U> class FluidDescriptor>
void computeTurbulentDiffusivityFromScalar(
    T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid, MultiScalarField3D<T> &D0,
    MultiScalarField3D<T> &diffusivity, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&fluid);
    args.push_back(&D0);
    args.push_back(&diffusivity);
    applyProcessingFunctional(
        new TurbulentDiffusivityFromScalarFunctional3D<T, FluidDescriptor>(Pr_t), domain, args);
}

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivityFromScalar(
    T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid, MultiScalarField3D<T> &D0, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > diffusivity =
        generateMultiScalarField<T>(fluid, domain);
    computeTurbulentDiffusivityFromScalar(Pr_t, fluid, D0, *diffusivity, domain);
    return diffusivity;
}

template <typename T, template <typename U> class FluidDescriptor>
std::unique_ptr<MultiScalarField3D<T> > computeTurbulentDiffusivityFromScalar(
    T Pr_t, MultiBlockLattice3D<T, FluidDescriptor> &fluid, MultiScalarField3D<T> &D0)
{
    return computeTurbulentDiffusivityFromScalar(Pr_t, fluid, D0, fluid.getBoundingBox());
}

template <typename T>
T ArrheniusChemicalReactionCoupling3D<T>::R = 8.3144598;

template <typename T>
ArrheniusChemicalReactionCoupling3D<T>::ArrheniusChemicalReactionCoupling3D(
    plint nSpecies_, plint nReactions_, T dt_, std::vector<T> const &preExponentialFactors_,
    std::vector<T> const &activationEnergies_, std::vector<T> const &reactionEnthalpies_,
    std::vector<std::vector<int> > const &stoichiometryMatrix_,
    std::vector<std::vector<T> > const &reactionRateMatrix_) :
    nSpecies(nSpecies_),
    nReactions(nReactions_),
    dt(dt_),
    preExponentialFactors(preExponentialFactors_),
    activationEnergies(activationEnergies_),
    reactionEnthalpies(reactionEnthalpies_),
    stoichiometryMatrix(stoichiometryMatrix_),
    reactionRateMatrix(reactionRateMatrix_)
{ }

template <typename T>
void ArrheniusChemicalReactionCoupling3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() >= 5);

    plint iBlock = 0;

    ScalarField3D<T> *theta = dynamic_cast<ScalarField3D<T> *>(blocks[iBlock]);
    PLB_ASSERT(theta);
    Dot3D ofsTheta(0, 0, 0);
    iBlock++;

    ScalarField3D<T> *thermalQ = dynamic_cast<ScalarField3D<T> *>(blocks[iBlock]);
    PLB_ASSERT(thermalQ);
    Dot3D ofsTQ = computeRelativeDisplacement(*theta, *thermalQ);
    iBlock++;

    std::vector<ScalarField3D<T> *> C(nSpecies);
    std::vector<Dot3D> ofsC(nSpecies);
    for (plint iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        C[iSpecies] = dynamic_cast<ScalarField3D<T> *>(blocks[iBlock]);
        PLB_ASSERT(C[iSpecies]);
        ofsC[iSpecies] = computeRelativeDisplacement(*theta, *C[iSpecies]);
        iBlock++;
    }

    std::vector<ScalarField3D<T> *> massQ(nSpecies);
    std::vector<Dot3D> ofsMQ(nSpecies);
    for (plint iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        massQ[iSpecies] = dynamic_cast<ScalarField3D<T> *>(blocks[iBlock]);
        PLB_ASSERT(massQ[iSpecies]);
        ofsMQ[iSpecies] = computeRelativeDisplacement(*theta, *massQ[iSpecies]);
        iBlock++;
    }

    ScalarField3D<T> *VHC = dynamic_cast<ScalarField3D<T> *>(blocks[iBlock]);
    PLB_ASSERT(VHC);
    Dot3D ofsVHC = computeRelativeDisplacement(*theta, *VHC);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // Reaction rates.

                std::vector<T> r(nReactions);
                for (plint iReaction = 0; iReaction < nReactions; iReaction++) {
                    r[iReaction] =
                        preExponentialFactors[iReaction]
                        * std::exp(
                            -activationEnergies[iReaction]
                            / (R * theta->get(iX + ofsTheta.x, iY + ofsTheta.y, iZ + ofsTheta.z)));
                    for (plint iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                        r[iReaction] *= std::pow(
                            C[iSpecies]->get(
                                iX + ofsC[iSpecies].x, iY + ofsC[iSpecies].y,
                                iZ + ofsC[iSpecies].z),
                            reactionRateMatrix[iReaction][iSpecies]);
                    }
                }

                // Source for the heat equation.

                T thetaQ = (T)0;
                for (plint iReaction = 0; iReaction < nReactions; iReaction++) {
                    thetaQ += r[iReaction] * reactionEnthalpies[iReaction];
                }
                thermalQ->get(iX + ofsTQ.x, iY + ofsTQ.y, iZ + ofsTQ.z) =
                    (thetaQ / VHC->get(iX + ofsVHC.x, iY + ofsVHC.y, iZ + ofsVHC.z)) * dt;

                // Sources for the concentration equations.

                for (plint iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
                    T CQ = (T)0;
                    for (plint iReaction = 0; iReaction < nReactions; iReaction++) {
                        CQ += stoichiometryMatrix[iSpecies][iReaction] * r[iReaction];
                    }
                    massQ[iSpecies]->get(
                        iX + ofsMQ[iSpecies].x, iY + ofsMQ[iSpecies].y, iZ + ofsMQ[iSpecies].z) =
                        CQ * dt;
                }
            }
        }
    }
}

template <typename T>
ArrheniusChemicalReactionCoupling3D<T> *ArrheniusChemicalReactionCoupling3D<T>::clone() const
{
    return new ArrheniusChemicalReactionCoupling3D<T>(*this);
}

template <typename T>
void ArrheniusChemicalReactionCoupling3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    PLB_ASSERT(modified.size() >= 5);

    plint iBlock = 0;

    modified[iBlock] = modif::nothing;  // Temperature.
    iBlock++;
    modified[iBlock] = modif::staticVariables;  // Temperature source.
    iBlock++;
    for (plint iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        modified[iBlock] = modif::nothing;  // Mass concentrations.
        iBlock++;
    }
    for (plint iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        modified[iBlock] = modif::staticVariables;  // Mass concentration sources.
        iBlock++;
    }
    modified[iBlock] = modif::nothing;  // Volumetric heat capacity.
}

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_3D_HH
