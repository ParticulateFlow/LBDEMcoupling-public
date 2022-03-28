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

#ifndef FREE_SURFACE_INITIALIZER_3D_HH
#define FREE_SURFACE_INITIALIZER_3D_HH

#include <cstdlib>

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/globalDefs.h"
#include "multiPhysics/freeSurfaceInitializer3D.h"
#include "sitmo/prng_engine.hpp"

namespace plb {

/* *************** Class DefaultInitializeFreeSurface3D *******************************************
 */

template <typename T, template <typename U> class Descriptor>
void DefaultInitializeFreeSurface3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    typedef Descriptor<T> D;

    // In the following, spot the interface cells and tag them.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (isEmpty(param.flag(iX, iY, iZ))) {
                    for (plint iPop = 1; iPop < D::q; iPop++) {
                        plint nextX = iX + D::c[iPop][0];
                        plint nextY = iY + D::c[iPop][1];
                        plint nextZ = iZ + D::c[iPop][2];
                        // Note: there is no conflict of concurrent read/write on param.flag,
                        // because read tests is a cell is fluid, and write only converts empty
                        // cells to interface.
                        if (isFullWet(param.flag(nextX, nextY, nextZ))) {
                            param.flag(iX, iY, iZ) = interface;
                        }
                    }
                }
            }
        }
    }
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T rho = 0.;
                if (useRhoIni) {
                    rho = rhoIni;
                } else {
                    rho = param.getDensity(iX, iY, iZ);
                }
                Array<T, 3> j((T)0., (T)0., (T)0.);
                if (!useZeroMomentum) {
                    j = param.getMomentum(iX, iY, iZ);
                }

                if (initializeCell) {
                    T rhoBar = Descriptor<T>::rhoBar(rho);
                    T jSqr = normSqr(j);
                    T thetaBar = 0.;
                    param.cell(iX, iY, iZ)
                        .getDynamics()
                        .computeEquilibria(
                            param.cell(iX, iY, iZ).getRawPopulations(), rhoBar, j, jSqr, thetaBar);
                }

                param.setDensity(iX, iY, iZ, rho);
                param.setMomentum(iX, iY, iZ, j);
                param.setForce(iX, iY, iZ, force);
                param.outsideDensity(iX, iY, iZ) = rhoIni;
                switch (param.flag(iX, iY, iZ)) {
                case fluid:
                case protect:
                case temporarilyProtect:
                    param.attributeDynamics(iX, iY, iZ, dynamicsTemplate->clone());
                    param.mass(iX, iY, iZ) = rho;
                    param.volumeFraction(iX, iY, iZ) = (T)1.;
                    break;
                case interface:
                    param.attributeDynamics(iX, iY, iZ, dynamicsTemplate->clone());
                    param.mass(iX, iY, iZ) = 0.5 * rho;
                    param.volumeFraction(iX, iY, iZ) = (T)0.5;
                    break;
                case empty:
                case protectEmpty:
                    param.attributeDynamics(iX, iY, iZ, emptyNodeDynamicsTemplate->clone());
                    param.setForce(iX, iY, iZ, Array<T, D::ExternalField::sizeOfForce>::zero());
                    param.mass(iX, iY, iZ) = (T)0.;
                    param.volumeFraction(iX, iY, iZ) = (T)0.;
                    break;
                case wall:
                    param.attributeDynamics(iX, iY, iZ, new BounceBack<T, Descriptor>(rhoIni));
                    param.mass(iX, iY, iZ) = (T)0.;
                    param.volumeFraction(iX, iY, iZ) = (T)0.;
                    break;
                case slipWall:
                    // If no reflection plane is given, SpecularReflection behaves like BounceBack.
                    param.attributeDynamics(
                        iX, iY, iZ, new SpecularReflection<T, Descriptor>(rhoIni));
                    param.mass(iX, iY, iZ) = (T)0.;
                    param.volumeFraction(iX, iY, iZ) = (T)0.;
                    break;
                default:
                    // Invalid free-surface flag.
                    PLB_ASSERT(false);
                }
            }
        }
    }
}

/* *************** Class PartiallyDefaultInitializeFreeSurface3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void PartiallyDefaultInitializeFreeSurface3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    typedef Descriptor<T> D;

    // In the following, spot the interface cells and tag them.
    // This time set the volume fraction to 0.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (isEmpty(param.flag(iX, iY, iZ))) {
                    for (plint iPop = 1; iPop < D::q; iPop++) {
                        plint nextX = iX + D::c[iPop][0];
                        plint nextY = iY + D::c[iPop][1];
                        plint nextZ = iZ + D::c[iPop][2];
                        // Note: there is no conflict of concurrent read/write on param.flag,
                        // because read tests is a cell is fluid, and write only converts empty
                        // cells to interface.
                        if (isFullWet(param.flag(nextX, nextY, nextZ))) {
                            param.flag(iX, iY, iZ) = interface;
                            param.volumeFraction(iX, iY, iZ) = (T)0;
                        }
                    }
                }
            }
        }
    }
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T rho = 0.;
                if (useRhoIni) {
                    rho = rhoIni;
                } else {
                    rho = param.getDensity(iX, iY, iZ);
                }
                Array<T, 3> j((T)0., (T)0., (T)0.);
                if (!useZeroMomentum) {
                    j = param.getMomentum(iX, iY, iZ);
                }

                if (initializeCell) {
                    T rhoBar = Descriptor<T>::rhoBar(rho);
                    T jSqr = normSqr(j);
                    T thetaBar = 0.;
                    param.cell(iX, iY, iZ)
                        .getDynamics()
                        .computeEquilibria(
                            param.cell(iX, iY, iZ).getRawPopulations(), rhoBar, j, jSqr, thetaBar);
                }

                param.setDensity(iX, iY, iZ, rho);
                param.setMomentum(iX, iY, iZ, j);
                param.setForce(iX, iY, iZ, force);
                param.outsideDensity(iX, iY, iZ) = rhoIni;
                switch (param.flag(iX, iY, iZ)) {
                case fluid:
                case protect:
                case temporarilyProtect:
                    param.attributeDynamics(iX, iY, iZ, dynamicsTemplate->clone());
                    param.mass(iX, iY, iZ) = rho;
                    break;
                case interface:
                    param.attributeDynamics(iX, iY, iZ, dynamicsTemplate->clone());
                    param.mass(iX, iY, iZ) = rho * param.volumeFraction(iX, iY, iZ);
                    break;
                case empty:
                case protectEmpty:
                    param.attributeDynamics(iX, iY, iZ, emptyNodeDynamicsTemplate->clone());
                    param.setForce(iX, iY, iZ, Array<T, D::ExternalField::sizeOfForce>::zero());
                    param.mass(iX, iY, iZ) = (T)0.;
                    break;
                case wall:
                    param.attributeDynamics(iX, iY, iZ, new BounceBack<T, Descriptor>(rhoIni));
                    param.mass(iX, iY, iZ) = (T)0.;
                    break;
                case slipWall:
                    // If no reflection plane is given, SpecularReflection behaves like BounceBack.
                    param.attributeDynamics(
                        iX, iY, iZ, new SpecularReflection<T, Descriptor>(rhoIni));
                    param.mass(iX, iY, iZ) = (T)0.;
                    break;
                default:
                    // Invalid free-surface flag.
                    PLB_ASSERT(false);
                }
            }
        }
    }
}

/* *************** Class ConstantIniVelocityFreeSurface3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void ConstantIniVelocityFreeSurface3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    T rho = rhoIni;
    Array<T, 3> j(velocity * rho);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (isWet(param.flag(iX, iY, iZ))) {
                    iniCellAtEquilibrium(param.cell(iX, iY, iZ), rho, velocity);
                    param.setMomentum(iX, iY, iZ, j);
                }
            }
        }
    }
}

/* *************** Class InletConstVolumeFraction3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void InletConstVolumeFraction3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                param.mass(iX, iY, iZ) = param.getDensity(iX, iY, iZ) * volumeFraction;
            }
        }
    }
}

/* *************** Class MaskedInletConstVolumeFraction3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void MaskedInletConstVolumeFraction3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 3);
    ScalarField3D<T> *rhoBar = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[0]);
    PLB_ASSERT(rhoBar);
    ScalarField3D<T> *mass = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[1]);
    PLB_ASSERT(mass);
    ScalarField3D<int> *mask = dynamic_cast<ScalarField3D<int> *>(atomicBlocks[2]);
    PLB_ASSERT(mask);

    Dot3D ofsMass = computeRelativeDisplacement(*rhoBar, *mass);
    Dot3D ofsMask = computeRelativeDisplacement(*rhoBar, *mask);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (mask->get(iX + ofsMask.x, iY + ofsMask.y, iZ + ofsMask.z) == whichFlag) {
                    mass->get(iX + ofsMass.x, iY + ofsMass.y, iZ + ofsMass.z) =
                        volumeFraction * Descriptor<T>::fullRho(rhoBar->get(iX, iY, iZ));
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void maskedInletConstVolumeFraction3D(
    MultiScalarField3D<T> &rhoBar, MultiScalarField3D<T> &mass, MultiScalarField3D<int> &mask,
    T volumeFraction, int whichFlag, Box3D domain)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&rhoBar);
    args.push_back(&mass);
    args.push_back(&mask);
    applyProcessingFunctional(
        new MaskedInletConstVolumeFraction3D<T, Descriptor>(volumeFraction, whichFlag), domain,
        args);
}

/* *************** Class N_MaskedInletConstVolumeFraction3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void N_MaskedInletConstVolumeFraction3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 3);
    ScalarField3D<T> *rhoBar = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[0]);
    PLB_ASSERT(rhoBar);
    ScalarField3D<T> *mass = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[1]);
    PLB_ASSERT(mass);
    NTensorField3D<int> *mask = dynamic_cast<NTensorField3D<int> *>(atomicBlocks[2]);
    PLB_ASSERT(mask);
    PLB_ASSERT(mask->getNdim() == 1);

    Dot3D ofsMass = computeRelativeDisplacement(*rhoBar, *mass);
    Dot3D ofsMask = computeRelativeDisplacement(*rhoBar, *mask);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (*mask->get(iX + ofsMask.x, iY + ofsMask.y, iZ + ofsMask.z)) {
                    mass->get(iX + ofsMass.x, iY + ofsMass.y, iZ + ofsMass.z) =
                        volumeFraction * Descriptor<T>::fullRho(rhoBar->get(iX, iY, iZ));
                }
            }
        }
    }
}

/* *************** Class OutletMaximumVolumeFraction3D *******************************************
 */

template <typename T, template <typename U> class Descriptor>
void OutletMaximumVolumeFraction3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T maximumMass = param.getDensity(iX, iY, iZ) * volumeFraction;
                if (param.mass(iX, iY, iZ) > maximumMass) {
                    param.mass(iX, iY, iZ) = maximumMass;
                }
            }
        }
    }
}

/* *************** Class N_MaskedOutletMaximumVolumeFraction3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void N_MaskedOutletMaximumVolumeFraction3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 3);
    ScalarField3D<T> *rhoBar = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[0]);
    PLB_ASSERT(rhoBar);
    ScalarField3D<T> *mass = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[1]);
    PLB_ASSERT(mass);
    NTensorField3D<int> *mask = dynamic_cast<NTensorField3D<int> *>(atomicBlocks[2]);
    PLB_ASSERT(mask);
    PLB_ASSERT(mask->getNdim() == 1);

    Dot3D ofsMass = computeRelativeDisplacement(*rhoBar, *mass);
    Dot3D ofsMask = computeRelativeDisplacement(*rhoBar, *mask);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (*mask->get(iX + ofsMask.x, iY + ofsMask.y, iZ + ofsMask.z)) {
                    T maximumMass =
                        volumeFraction * Descriptor<T>::fullRho(rhoBar->get(iX, iY, iZ));
                    if (mass->get(iX + ofsMass.x, iY + ofsMass.y, iZ + ofsMass.z) > maximumMass) {
                        mass->get(iX + ofsMass.x, iY + ofsMass.y, iZ + ofsMass.z) = maximumMass;
                    }
                }
            }
        }
    }
}

/* *************** Class OutletVolumeFractionInRange3D *******************************************
 */

template <typename T, template <typename U> class Descriptor>
void OutletVolumeFractionInRange3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T maximumMass = param.getDensity(iX, iY, iZ) * maxFraction;
                T minimumMass = param.getDensity(iX, iY, iZ) * minFraction;
                if (param.mass(iX, iY, iZ) > maximumMass) {
                    param.mass(iX, iY, iZ) = maximumMass;
                } else if (param.mass(iX, iY, iZ) < minimumMass) {
                    param.mass(iX, iY, iZ) = minimumMass;
                }
            }
        }
    }
}

/* *************** Class OutletMaximumVolumeFraction2_3D *******************************************
 */

template <typename T, template <typename U> class Descriptor>
void OutletMaximumVolumeFraction2_3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T maximumMass = param.getDensity(iX, iY, iZ) * volumeFraction;
                if (param.mass(iX, iY, iZ) > maximumMass) {
                    param.mass(iX, iY, iZ) = maximumMass;
                }

                Cell<T, Descriptor> &cell = param.cell(iX, iY, iZ);
                T oldRhoBar;
                Array<T, 3> oldJ;
                momentTemplates<T, Descriptor>::get_rhoBar_j(cell, oldRhoBar, oldJ);
                T oldJsqr = normSqr(oldJ);
                T rhoBar = Descriptor<T>::rhoBar(param.getDensity(iX, iY, iZ));
                for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    T oldEq = cell.getDynamics().computeEquilibrium(iPop, oldRhoBar, oldJ, oldJsqr);
                    T newEq = cell.getDynamics().computeEquilibrium(iPop, rhoBar, oldJ, oldJsqr);
                    cell[iPop] += newEq - oldEq;
                }
            }
        }
    }
}

/* *************** Class NoSlipMaximumVolumeFraction3D *******************************************
 */

template <typename T, template <typename U> class Descriptor>
void NoSlipMaximumVolumeFraction3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T maximumMass = param.getDensity(iX, iY, iZ) * volumeFraction;
                if (param.mass(iX, iY, iZ) > maximumMass) {
                    param.setDensity(iX, iY, iZ, param.mass(iX, iY, iZ) / volumeFraction);
                }
            }
        }
    }
}

/* *************** Class PunchSphere3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void PunchSphere3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    typedef Descriptor<T> D;

    Dot3D offset = param.absOffset();
    Array<T, 3> localCenter(center - Array<T, 3>(offset.x, offset.y, offset.z));

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (normSqr(Array<T, 3>(iX, iY, iZ) - localCenter) < radius * radius) {
                    bool isBoundary = false;
                    for (plint iPop = 0; iPop < D::q; ++iPop) {
                        plint nextX = iX + D::c[iPop][0];
                        plint nextY = iY + D::c[iPop][1];
                        plint nextZ = iZ + D::c[iPop][2];
                        if (normSqr((Array<T, 3>(nextX, nextY, nextZ) - localCenter))
                            >= radius * radius) {
                            isBoundary = true;
                        }
                    }
                    if (isBoundary) {
                        param.flag(iX, iY, iZ) = interface;
                        param.volumeFraction(iX, iY, iZ) = (T)0.5;
                        param.mass(iX, iY, iZ) = 0.5 * rho0;
                        param.setDensity(iX, iY, iZ, rho0);
                        param.outsideDensity(iX, iY, iZ) = rho0;
                    } else {
                        param.flag(iX, iY, iZ) = empty;
                        param.attributeDynamics(iX, iY, iZ, emptyNodeDynamicsTemplate->clone());
                        param.mass(iX, iY, iZ) = T();
                        param.volumeFraction(iX, iY, iZ) = T();
                        param.setDensity(iX, iY, iZ, rho0);
                        param.setForce(iX, iY, iZ, Array<T, D::ExternalField::sizeOfForce>::zero());
                        param.setMomentum(iX, iY, iZ, Array<T, 3>(T(), T(), T()));
                        param.outsideDensity(iX, iY, iZ) = rho0;
                    }
                }
            }
        }
    }
}

/* *************** Class AnalyticalPunchSphere3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void AnalyticalPunchSphere3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    Dot3D offset = param.absOffset();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint globalX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint globalY = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint globalZ = iZ + offset.z;

                if (param.flag(iX, iY, iZ) == empty || isAnyWall(param.flag(iX, iY, iZ))) {
                    continue;
                }

                int nextFlag = 0;
                T nextVolumeFraction = 0.;
                subDomainVolumeFraction(globalX, globalY, globalZ, nextFlag, nextVolumeFraction);

                if (nextFlag == fluid) {
                    // Do nothing.
                } else if (nextFlag == empty) {
                    param.flag(iX, iY, iZ) = empty;
                    param.attributeDynamics(iX, iY, iZ, emptyNodeDynamicsTemplate->clone());
                    param.mass(iX, iY, iZ) = T();
                    param.volumeFraction(iX, iY, iZ) = T();
                    param.setDensity(iX, iY, iZ, rho0);
                    // param.setForce(iX,iY,iZ,
                    // Array<T,Descriptor<T>::ExternalField::sizeOfForce>::zero());
                    param.setMomentum(iX, iY, iZ, Array<T, 3>(T(), T(), T()));
                    param.outsideDensity(iX, iY, iZ) = rho0;
                } else {
                    param.flag(iX, iY, iZ) = interface;
                    param.volumeFraction(iX, iY, iZ) = nextVolumeFraction;
                    param.mass(iX, iY, iZ) = nextVolumeFraction * rho0;
                    param.setDensity(iX, iY, iZ, rho0);
                    param.outsideDensity(iX, iY, iZ) = rho0;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void AnalyticalPunchSphere3D<T, Descriptor>::subDomainVolumeFraction(
    plint globalX, plint globalY, plint globalZ, int &flag, T &volumeFraction)
{
    plint numInside = 0;
    plint numOutside = 0;

    sitmo::prng_engine eng;
    T xi = (T)globalX - 0.5;
    T yi = (T)globalY - 0.5;
    T zi = (T)globalZ - 0.5;
    for (plint xSub = 0; xSub < subDivision; ++xSub) {
        T xPos = xi + (T)eng() / (T)sitmo::prng_engine::max();
        for (plint ySub = 0; ySub < subDivision; ++ySub) {
            T yPos = yi + (T)eng() / (T)sitmo::prng_engine::max();
            for (plint zSub = 0; zSub < subDivision; ++zSub) {
                T zPos = zi + (T)eng() / (T)sitmo::prng_engine::max();
                if (isInsideSphere(xPos, yPos, zPos)) {
                    ++numInside;
                } else {
                    ++numOutside;
                }
            }
        }
    }

    if (numInside == 0) {
        flag = freeSurfaceFlag::fluid;
        volumeFraction = (T)1;
    } else if (numOutside == 0) {
        flag = freeSurfaceFlag::empty;
        volumeFraction = (T)0;
    } else {
        flag = freeSurfaceFlag::interface;
        volumeFraction = (T)numOutside / ((T)numInside + (T)numOutside);
    }
}

/* *************** Class CalculateAverageSphereDensity3D *******************************************
 */

template <typename T, template <typename U> class Descriptor>
void CalculateAverageSphereDensity3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    Dot3D offset = param.absOffset();
    Array<T, 3> localCenter(center - Array<T, 3>(offset.x, offset.y, offset.z));
    BlockStatistics &statistics = this->getStatistics();

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (normSqr(Array<T, 3>(iX, iY, iZ) - localCenter) < radius * radius) {
                    statistics.gatherAverage(averageDensityId, param.getDensity(iX, iY, iZ));
                    statistics.incrementStats();
                }
            }
        }
    }
}

/* *************** Class AnalyticalIniVolumeFraction3D *******************************************
 */

template <typename T, class InsideFunction>
void AnalyticalIniVolumeFraction3D<T, InsideFunction>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 2);
    ScalarField3D<T> *volumeFraction = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[0]);
    PLB_ASSERT(volumeFraction);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int> *>(atomicBlocks[1]);
    PLB_ASSERT(flag);

    Dot3D offset = computeRelativeDisplacement(*volumeFraction, *flag);
    Dot3D absOfs = volumeFraction->getLocation();

    // In the following, spot the interface cells and tag them.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                int nextFlag = 0;
                T nextVolumeFraction = 0.;
                subDomainVolumeFraction(
                    iX + absOfs.x, iY + absOfs.y, iZ + absOfs.z, nextFlag, nextVolumeFraction);
                volumeFraction->get(iX, iY, iZ) = nextVolumeFraction;
                flag->get(iX + offset.x, iY + offset.y, iZ + offset.z) = nextFlag;
            }
        }
    }
}

template <typename T, class InsideFunction>
void AnalyticalIniVolumeFraction3D<T, InsideFunction>::subDomainVolumeFraction(
    plint iX, plint iY, plint iZ, int &flag, T &volumeFraction)
{
    plint numInside = 0;
    plint numOutside = 0;

    sitmo::prng_engine eng;
    T xi = (T)iX - 0.5;
    T yi = (T)iY - 0.5;
    T zi = (T)iZ - 0.5;
    for (plint xSub = 0; xSub < subDivision; ++xSub) {
        T xPos = xi + (T)eng() / (T)sitmo::prng_engine::max();
        for (plint ySub = 0; ySub < subDivision; ++ySub) {
            T yPos = yi + (T)eng() / (T)sitmo::prng_engine::max();
            for (plint zSub = 0; zSub < subDivision; ++zSub) {
                T zPos = zi + (T)eng() / (T)sitmo::prng_engine::max();
                if (insideFunction(xPos, yPos, zPos)) {
                    ++numInside;
                } else {
                    ++numOutside;
                }
            }
        }
    }

    if (numInside == 0) {
        flag = freeSurfaceFlag::empty;
        volumeFraction = (T)0;
    } else if (numOutside == 0) {
        flag = freeSurfaceFlag::fluid;
        volumeFraction = (T)1;
    } else {
        flag = freeSurfaceFlag::interface;
        volumeFraction = (T)numInside / ((T)numInside + (T)numOutside);
    }
}

template <typename T, class InsideFunction>
void analyticalIniVolumeFraction(
    MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flagStatus,
    InsideFunction const &insideFunction, Box3D domain, plint subDivision)
{
    std::vector<MultiBlock3D *> args;
    args.push_back(&volumeFraction);
    args.push_back(&flagStatus);
    applyProcessingFunctional(
        new AnalyticalIniVolumeFraction3D<T, InsideFunction>(insideFunction, subDivision), domain,
        args);
}

template <typename T, class InsideFunction>
void analyticalIniVolumeFraction(
    MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flagStatus,
    InsideFunction const &insideFunction, plint subDivision)
{
    analyticalIniVolumeFraction(
        volumeFraction, flagStatus, insideFunction, volumeFraction.getBoundingBox(), subDivision);
}

}  // namespace plb

#endif  // FREE_SURFACE_INITIALIZER_3D_HH
