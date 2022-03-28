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

#ifndef FREE_SURFACE_BOUNDARY_CONDITION_3D_HH
#define FREE_SURFACE_BOUNDARY_CONDITION_3D_HH

#include <cmath>
#include <iostream>

#include "multiPhysics/freeSurfaceBoundaryCondition3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
FreeSurfaceFadingArea3D<T, Descriptor>::FreeSurfaceFadingArea3D(T factor_) : factor(factor_)
{ }

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceFadingArea3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice)
{
    std::vector<T> decomposedVariables;

    enum {
        momentumStoredOffset = Descriptor<T>::ExternalField::momentumBeginsAt,
        densityStoredOffset = Descriptor<T>::ExternalField::densityBeginsAt
    };

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = lattice.get(iX, iY, iZ);
                plint order = 0;
                cell.getDynamics().decompose(cell, decomposedVariables, order);

                T density = Descriptor<T>::fullRho(decomposedVariables[0]);
                if (density > T(0))
                    density *= factor;
                decomposedVariables[0] = Descriptor<T>::rhoBar(density);
                cell.getDynamics().recompose(cell, decomposedVariables, order);

                *cell.getExternal(densityStoredOffset) = density;

                Array<T, Descriptor<T>::d> j;
                j.resetToZero();
                T rhoBar;
                momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);

                // TODO: What about mass, volumeFraction, flagStatus?
                j.to_cArray(cell.getExternal(momentumStoredOffset));
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceFadingArea3D<T, Descriptor> *FreeSurfaceFadingArea3D<T, Descriptor>::clone() const
{
    return new FreeSurfaceFadingArea3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void RemoveMass3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // param.attributeDynamics(iX,iY,iZ, new NoDynamics<T,Descriptor>((T)1.));
                param.setDensity(iX, iY, iZ, (T)1.);
                param.setMomentum(iX, iY, iZ, Array<T, 3>((T)0., (T)0., (T)0.));
                param.mass(iX, iY, iZ) = (T)0;
                param.volumeFraction(iX, iY, iZ) = (T)0;
                // param.flag(iX,iY,iZ) = empty;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
RemoveMass3D<T, Descriptor> *RemoveMass3D<T, Descriptor>::clone() const
{
    return new RemoveMass3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void PouringLiquid3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T iniRho = T(1);
                param.attributeDynamics(iX, iY, iZ, dynamicsTemplate->clone());
                iniCellAtEquilibrium(param.cell(iX, iY, iZ), iniRho, injectionVelocity);
                param.setDensity(iX, iY, iZ, iniRho);
                param.setMomentum(iX, iY, iZ, iniRho * injectionVelocity);
                param.mass(iX, iY, iZ) = iniRho;
                param.volumeFraction(iX, iY, iZ) = (T)1;
                param.flag(iX, iY, iZ) = fluid;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
PouringLiquid3D<T, Descriptor> *PouringLiquid3D<T, Descriptor>::clone() const
{
    return new PouringLiquid3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ShortenBounceBack3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    typedef Descriptor<T> D;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    Box3D extDomain = domain.enlarge(1);

    for (plint iX = extDomain.x0; iX <= extDomain.x1; ++iX) {
        bool xBoundary = iX == extDomain.x0 || iX == extDomain.x1;
        for (plint iY = extDomain.y0; iY <= extDomain.y1; ++iY) {
            bool yBoundary = xBoundary || iY == extDomain.y0 || iY == extDomain.y1;
            for (plint iZ = extDomain.z0; iZ <= extDomain.z1; ++iZ) {
                if (isWall(param.flag(iX, iY, iZ)))
                {  // Here we want only the no-slip wall cells (not the free-slip ones).
                    bool zBoundary = yBoundary || iZ == extDomain.z0 || iZ == extDomain.z1;
                    for (plint iNeighbor = 1; iNeighbor < D::q; ++iNeighbor) {
                        plint nextX = iX + D::c[iNeighbor][0];
                        plint nextY = iY + D::c[iNeighbor][1];
                        plint nextZ = iZ + D::c[iNeighbor][2];
                        if (!zBoundary || contained(nextX, nextY, nextZ, domain)) {
                            if (isWet(param.flag(nextX, nextY, nextZ))) {
                                plint opp = indexTemplates::opposite<D>(iNeighbor);
                                param.cell(nextX, nextY, nextZ)[iNeighbor] =
                                    param.cell(iX, iY, iZ)[opp];
                            }
                        }
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ShortenBounceBack3D<T, Descriptor> *ShortenBounceBack3D<T, Descriptor>::clone() const
{
    return new ShortenBounceBack3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceSpongeZone3D<T, Descriptor>::FreeSurfaceSpongeZone3D(
    plint nx_, plint ny_, plint nz_, Array<plint, 6> const &numSpongeCells_,
    Array<T, 6> const &translationParameters_, Array<T, 6> const &scaleParameters_,
    bool incompressibleModel_) :
    nx(nx_),
    ny(ny_),
    nz(nz_),
    numSpongeCells(numSpongeCells_),
    translationParameters(translationParameters_),
    scaleParameters(scaleParameters_),
    incompressibleModel(incompressibleModel_),
    useTanhSpongeFunction(true)
{
    for (int iZone = 0; iZone < 6; iZone++) {
        if (numSpongeCells[iZone] > 0) {
            PLB_ASSERT(translationParameters[iZone] > (T)0 && translationParameters[iZone] < (T)1);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceSpongeZone3D<T, Descriptor>::FreeSurfaceSpongeZone3D(
    plint nx_, plint ny_, plint nz_, Array<plint, 6> const &numSpongeCells_,
    bool incompressibleModel_) :
    nx(nx_),
    ny(ny_),
    nz(nz_),
    numSpongeCells(numSpongeCells_),
    incompressibleModel(incompressibleModel_),
    useTanhSpongeFunction(false)
{ }

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceSpongeZone3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 6);

    ScalarField3D<T> *volumeFraction = dynamic_cast<ScalarField3D<T> *>(blocks[0]);
    PLB_ASSERT(volumeFraction);
    ScalarField3D<T> *targetVolumeFraction = dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    PLB_ASSERT(targetVolumeFraction);
    ScalarField3D<T> *rhoBar = dynamic_cast<ScalarField3D<T> *>(blocks[2]);
    PLB_ASSERT(rhoBar);
    ScalarField3D<T> *mass = dynamic_cast<ScalarField3D<T> *>(blocks[3]);
    PLB_ASSERT(mass);
    TensorField3D<T, 3> *j = dynamic_cast<TensorField3D<T, 3> *>(blocks[4]);
    PLB_ASSERT(j);
    TensorField3D<T, 3> *targetVelocity = dynamic_cast<TensorField3D<T, 3> *>(blocks[5]);
    PLB_ASSERT(targetVelocity);

    Dot3D location = volumeFraction->getLocation();

    Dot3D ofsTVF = computeRelativeDisplacement(*volumeFraction, *targetVolumeFraction);
    Dot3D ofsRB = computeRelativeDisplacement(*volumeFraction, *rhoBar);
    Dot3D ofsM = computeRelativeDisplacement(*volumeFraction, *mass);
    Dot3D ofsJ = computeRelativeDisplacement(*volumeFraction, *j);
    Dot3D ofsTV = computeRelativeDisplacement(*volumeFraction, *targetVelocity);

    plint spongePositions[6] = {0};
    spongePositions[0] = numSpongeCells[0];
    spongePositions[1] = nx - 1 - numSpongeCells[1];
    spongePositions[2] = numSpongeCells[2];
    spongePositions[3] = ny - 1 - numSpongeCells[3];
    spongePositions[4] = numSpongeCells[4];
    spongePositions[5] = nz - 1 - numSpongeCells[5];

    plint spongeCenters[6] = {0};
    T invSigma[6];
    T invNumSpongeCells[6];
    if (useTanhSpongeFunction) {
        spongeCenters[0] = util::roundToInt(((T)1 - translationParameters[0]) * numSpongeCells[0]);
        spongeCenters[1] =
            util::roundToInt(spongePositions[1] + translationParameters[1] * numSpongeCells[1]);
        spongeCenters[2] = util::roundToInt(((T)1 - translationParameters[2]) * numSpongeCells[2]);
        spongeCenters[3] =
            util::roundToInt(spongePositions[3] + translationParameters[3] * numSpongeCells[3]);
        spongeCenters[4] = util::roundToInt(((T)1 - translationParameters[4]) * numSpongeCells[4]);
        spongeCenters[5] =
            util::roundToInt(spongePositions[5] + translationParameters[5] * numSpongeCells[5]);

        for (int iZone = 0; iZone < 6; iZone++) {
            invSigma[iZone] = (T)1 / (scaleParameters[iZone] * (T)numSpongeCells[iZone]);
        }
    } else {
        for (int iZone = 0; iZone < 6; iZone++) {
            invNumSpongeCells[iZone] = (T)1 / (T)numSpongeCells[iZone];
        }
    }
    T pi = std::acos((T)-1);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint z = iZ + location.z;

                T spongeFunction = (T)1;
                bool modify = false;
                if (numSpongeCells[0] > 0 && x <= spongePositions[0]) {
                    if (useTanhSpongeFunction) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(x - spongeCenters[0]) * invSigma[0]));
                    } else {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)x * invNumSpongeCells[0]));
                    }
                    modify = true;
                }
                if (numSpongeCells[1] > 0 && x >= spongePositions[1]) {
                    if (useTanhSpongeFunction) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(x - spongeCenters[1]) * invSigma[1]));
                    } else {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(x - spongePositions[1]) * invNumSpongeCells[1]));
                    }
                    modify = true;
                }

                if (numSpongeCells[2] > 0 && y <= spongePositions[2]) {
                    if (useTanhSpongeFunction) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(y - spongeCenters[2]) * invSigma[2]));
                    } else {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)y * invNumSpongeCells[2]));
                    }
                    modify = true;
                }
                if (numSpongeCells[3] > 0 && y >= spongePositions[3]) {
                    if (useTanhSpongeFunction) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(y - spongeCenters[3]) * invSigma[3]));
                    } else {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(y - spongePositions[3]) * invNumSpongeCells[3]));
                    }
                    modify = true;
                }

                if (numSpongeCells[4] > 0 && z <= spongePositions[4]) {
                    if (useTanhSpongeFunction) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(z - spongeCenters[4]) * invSigma[4]));
                    } else {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)z * invNumSpongeCells[4]));
                    }
                    modify = true;
                }
                if (numSpongeCells[5] > 0 && z >= spongePositions[5]) {
                    if (useTanhSpongeFunction) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(z - spongeCenters[5]) * invSigma[5]));
                    } else {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(z - spongePositions[5]) * invNumSpongeCells[5]));
                    }
                    modify = true;
                }

                if (modify) {
                    T rate = (T)1 - spongeFunction;

                    T currentVF = volumeFraction->get(iX, iY, iZ);
                    T targetVF =
                        targetVolumeFraction->get(iX + ofsTVF.x, iY + ofsTVF.y, iZ + ofsTVF.z);
                    T newVolumeFraction = currentVF + rate * (targetVF - currentVF);
                    volumeFraction->get(iX, iY, iZ) = newVolumeFraction;

                    T rho = Descriptor<T>::fullRho(
                        rhoBar->get(iX + ofsRB.x, iY + ofsRB.y, iZ + ofsRB.z));
                    mass->get(iX + ofsM.x, iY + ofsM.y, iZ + ofsM.z) = rho * newVolumeFraction;

                    Array<T, 3> const &targetV =
                        targetVelocity->get(iX + ofsTV.x, iY + ofsTV.y, iZ + ofsTV.z);
                    Array<T, 3> const &currentJ = j->get(iX + ofsJ.x, iY + ofsJ.y, iZ + ofsJ.z);
                    if (incompressibleModel) {
                        j->get(iX + ofsJ.x, iY + ofsJ.y, iZ + ofsJ.z) =
                            currentJ + rate * (targetV - currentJ);
                    } else {
                        j->get(iX + ofsJ.x, iY + ofsJ.y, iZ + ofsJ.z) =
                            currentJ + rate * (rho * targetV - currentJ);
                    }
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
FreeSurfaceSpongeZone3D<T, Descriptor> *FreeSurfaceSpongeZone3D<T, Descriptor>::clone() const
{
    return new FreeSurfaceSpongeZone3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceSpongeZone3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Volume fraction.
    modified[1] = modif::nothing;          // Target volume fraction.
    modified[2] = modif::nothing;          // rhoBar.
    modified[3] = modif::staticVariables;  // mass.
    modified[4] = modif::staticVariables;  // j.
    modified[5] = modif::nothing;          // Target j.
}

}  // namespace plb

#endif  // FREE_SURFACE_BOUNDARY_CONDITION_3D_HH
