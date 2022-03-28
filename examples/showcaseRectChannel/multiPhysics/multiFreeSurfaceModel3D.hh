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

#ifndef MULTI_FREE_SURFACE_MODEL_3D_HH
#define MULTI_FREE_SURFACE_MODEL_3D_HH

#include <cmath>

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiPhysics/freeSurfaceModel3D.h"
#include "multiPhysics/freeSurfaceTemplates.h"
#include "multiPhysics/multiFreeSurfaceModel3D.h"

namespace plb {

/* *************** Class MultiFreeSurfaceOneWayCoupling3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void MultiFreeSurfaceOneWayCoupling3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 20);

    using namespace freeSurfaceFlag;

    std::vector<AtomicBlock3D *> atomicBlocks1, atomicBlocks2;

    for (int i = 0; i < 10; i++)
        atomicBlocks1.push_back(atomicBlocks[i]);

    for (int i = 0; i < 10; i++)
        atomicBlocks2.push_back(atomicBlocks[10 + i]);

    FreeSurfaceProcessorParam3D<T, Descriptor> param1(atomicBlocks1);
    FreeSurfaceProcessorParam3D<T, Descriptor> param2(atomicBlocks2);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (param1.flag(iX, iY, iZ) == interface && isWet(param2.flag(iX, iY, iZ))) {
                    // Velocity coupling
                    T rhoBar1;
                    Array<T, 3> j1;
                    momentTemplates<T, Descriptor>::get_rhoBar_j(
                        param1.cell(iX, iY, iZ), rhoBar1, j1);
                    T density1 = Descriptor<T>::fullRho(rhoBar1);

                    T rhoBar2;
                    Array<T, 3> j2;
                    momentTemplates<T, Descriptor>::get_rhoBar_j(
                        param2.cell(iX, iY, iZ), rhoBar2, j2);
                    T density2 = Descriptor<T>::fullRho(rhoBar2);

                    T dynamicOmega1 =
                        param1.cell(iX, iY, iZ)
                            .getDynamics()
                            .getDynamicParameter(
                                dynamicParams::dynamicOmega,
                                param1.cell(iX, iY, iZ));  // In case of a Smagorinsky model.
                    T tau1 = 0.0;
                    if (!util::isZero(dynamicOmega1)) {
                        tau1 = (T)1 / dynamicOmega1;
                    } else {
                        tau1 = (T)1 / param1.cell(iX, iY, iZ).getDynamics().getOmega();
                    }

                    for (plint iD = 0; iD < 3; ++iD) {
                        // j1[iD] += interactionStrength*density1*(j2[iD]/density2 -
                        // j1[iD]/density1)*tau1;
                        j1[iD] += interactionStrength * rhoDefault1
                                  * (j2[iD] / density2 - j1[iD] / density1) * tau1;
                    }

                    param1.setMomentum(iX, iY, iZ, j1);
                }
            }
        }
    }
}

/* *************** Class MultiFreeSurfaceVelocityContinuityCoupling3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void MultiFreeSurfaceVelocityContinuityCoupling3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 20);

    using namespace freeSurfaceFlag;

    std::vector<AtomicBlock3D *> atomicBlocks1, atomicBlocks2;

    for (int i = 0; i < 10; i++)
        atomicBlocks1.push_back(atomicBlocks[i]);

    for (int i = 0; i < 10; i++)
        atomicBlocks2.push_back(atomicBlocks[10 + i]);

    FreeSurfaceProcessorParam3D<T, Descriptor> param1(atomicBlocks1);
    FreeSurfaceProcessorParam3D<T, Descriptor> param2(atomicBlocks2);

    T convertRho1to2 = rhoDefault2 / rhoDefault1;
    T convertRho2to1 = rhoDefault1 / rhoDefault2;

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 3> velocity((T)0.0, (T)0.0, (T)0.0);

                // Enforce continuity of velocity only when both fluid cells belong to the
                // interface.
                if (param1.flag(iX, iY, iZ) == interface && param2.flag(iX, iY, iZ) == interface) {
                    // if ((param2.flag(iX,iY,iZ)==interface && isWet(param1.flag(iX,iY,iZ))) ||
                    //     (param1.flag(iX,iY,iZ)==interface && isWet(param2.flag(iX,iY,iZ)))) {
                    T rhoBar1;
                    Array<T, 3> j1;
                    momentTemplates<T, Descriptor>::get_rhoBar_j(
                        param1.cell(iX, iY, iZ), rhoBar1, j1);
                    T density1 = Descriptor<T>::fullRho(rhoBar1);

                    T rhoBar2;
                    Array<T, 3> j2;
                    momentTemplates<T, Descriptor>::get_rhoBar_j(
                        param2.cell(iX, iY, iZ), rhoBar2, j2);
                    T density2 = Descriptor<T>::fullRho(rhoBar2);

                    velocity = 0.5 * (j1 / density1 + j2 / density2);
                }

                if (param1.flag(iX, iY, iZ) == interface && isWet(param2.flag(iX, iY, iZ))) {
                    // Compute the average pressure from fluid cells of both fluids.
                    T averageDensity = T();
                    plint numNeighbors = 0;
                    for (plint i = -1; i < 2; i++) {
                        for (plint j = -1; j < 2; j++) {
                            for (plint k = -1; k < 2; k++) {
                                if (isFullWet(param1.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity += param1.getDensity(iX + i, iY + j, iZ + k);
                                }
                                if (isFullWet(param2.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity +=
                                        convertRho2to1 * param2.getDensity(iX + i, iY + j, iZ + k);
                                }
                            }
                        }
                    }
                    if (numNeighbors == 0) {
                        averageDensity = rhoDefault1;
                    } else {
                        averageDensity /= (T)numNeighbors;
                    }
                    param1.setDensity(iX, iY, iZ, averageDensity);

                    param1.volumeFraction(iX, iY, iZ) = param1.mass(iX, iY, iZ) / averageDensity;
                    if (param2.flag(iX, iY, iZ) == interface) {
                        param1.setMomentum(iX, iY, iZ, param1.getDensity(iX, iY, iZ) * velocity);
                    }
                }

                if (param2.flag(iX, iY, iZ) == interface && isWet(param1.flag(iX, iY, iZ))) {
                    // Compute the average pressure from fluid cells of both fluids.
                    T averageDensity = T();
                    plint numNeighbors = 0;
                    for (plint i = -1; i < 2; i++) {
                        for (plint j = -1; j < 2; j++) {
                            for (plint k = -1; k < 2; k++) {
                                if (isFullWet(param1.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity +=
                                        convertRho1to2 * param1.getDensity(iX + i, iY + j, iZ + k);
                                }
                                if (isFullWet(param2.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity += param2.getDensity(iX + i, iY + j, iZ + k);
                                }
                            }
                        }
                    }
                    if (numNeighbors == 0) {
                        averageDensity = rhoDefault2;
                    } else {
                        averageDensity /= (T)numNeighbors;
                    }
                    param2.setDensity(iX, iY, iZ, averageDensity);

                    param2.volumeFraction(iX, iY, iZ) = param2.mass(iX, iY, iZ) / averageDensity;
                    if (param1.flag(iX, iY, iZ) == interface) {
                        param2.setMomentum(iX, iY, iZ, param2.getDensity(iX, iY, iZ) * velocity);
                    }
                }
            }
        }
    }
}

/* *************** Class MultiFreeSurfaceRepellingForceCoupling3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
T MultiFreeSurfaceRepellingForceCoupling3D<T, Descriptor>::deltaFunction(T r, T h)
{
    PLB_ASSERT(util::greaterEqual(r, (T)0));
    PLB_ASSERT(util::greaterThan(h, (T)0));

    static T pi = std::acos((T)-1);

    T delta = 0.0;
    if (util::lessEqual(r, h)) {
        delta = 0.5 * (1.0 + std::cos(pi * r / h));
    }

    return delta;
}

template <typename T, template <typename U> class Descriptor>
void MultiFreeSurfaceRepellingForceCoupling3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 20);

    using namespace freeSurfaceFlag;

    std::vector<AtomicBlock3D *> atomicBlocks1, atomicBlocks2;

    for (int i = 0; i < 10; i++)
        atomicBlocks1.push_back(atomicBlocks[i]);

    for (int i = 0; i < 10; i++)
        atomicBlocks2.push_back(atomicBlocks[10 + i]);

    FreeSurfaceProcessorParam3D<T, Descriptor> param1(atomicBlocks1);
    FreeSurfaceProcessorParam3D<T, Descriptor> param2(atomicBlocks2);

    T convertRho1to2 = rhoDefault2 / rhoDefault1;
    T convertRho2to1 = rhoDefault1 / rhoDefault2;

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // Implement the repelling force in terms of a momentum correction. Remember that
                // fluids can have strong density differences. In order to have the units of a
                // momentum, this correction term must be multiplied by the density.
                if (isWet(param1.flag(iX, iY, iZ)) && isWet(param2.flag(iX, iY, iZ))) {
                    // Add a force to both fluids. Remember the envelope is 3.
                    Array<T, 3> force1((T)0, (T)0, (T)0);
                    for (plint i = -2; i < 3; i++) {
                        plint I = iX + i;
                        for (plint j = -2; j < 3; j++) {
                            plint J = iY + j;
                            for (plint k = -2; k < 3; k++) {
                                plint K = iZ + k;
                                if (isWet(param1.flag(I, J, K)) && isWet(param2.flag(I, J, K))) {
                                    T VFdiff = param1.volumeFraction(I, J, K)
                                               + param2.volumeFraction(I, J, K) - 1.0;
                                    if (VFdiff > (T)0.0) {
                                        Array<T, 3> outwardNormalToFluid1 =
                                            param1.getNormal(I, J, K) - param2.getNormal(I, J, K);
                                        outwardNormalToFluid1 /= norm(outwardNormalToFluid1);
                                        T r = std::sqrt((T)i * i + (T)j * j + (T)k * k);
                                        T delta = deltaFunction(
                                            r, (T)3);  // h is 3 because envelope is 3.
                                        force1 += -VFdiff * delta * outwardNormalToFluid1;
                                    }
                                }
                            }
                        }
                    }
                    force1 *= interactionStrength;

                    T dynamicOmega1 =
                        param1.cell(iX, iY, iZ)
                            .getDynamics()
                            .getDynamicParameter(
                                dynamicParams::dynamicOmega,
                                param1.cell(iX, iY, iZ));  // In case of a Smagorinsky model.
                    T tau1 = 0.0;
                    if (!util::isZero(dynamicOmega1)) {
                        tau1 = (T)1 / dynamicOmega1;
                    } else {
                        tau1 = (T)1 / param1.cell(iX, iY, iZ).getDynamics().getOmega();
                    }

                    T dynamicOmega2 =
                        param2.cell(iX, iY, iZ)
                            .getDynamics()
                            .getDynamicParameter(
                                dynamicParams::dynamicOmega,
                                param2.cell(iX, iY, iZ));  // In case of a Smagorinsky model.
                    T tau2 = 0.0;
                    if (!util::isZero(dynamicOmega2)) {
                        tau2 = (T)1 / dynamicOmega2;
                    } else {
                        tau2 = (T)1 / param2.cell(iX, iY, iZ).getDynamics().getOmega();
                    }

                    // force1 is an acceleration.
                    param1.setMomentum(
                        iX, iY, iZ, param1.getMomentum(iX, iY, iZ) + rhoDefault1 * tau1 * force1);
                    param2.setMomentum(
                        iX, iY, iZ,
                        param2.getMomentum(iX, iY, iZ) + rhoDefault2 * tau2 * (-force1));

                    /*
                    // force1 is a force.
                    param1.setMomentum(iX,iY,iZ, param1.getMomentum(iX,iY,iZ) + tau1*force1);
                    param2.setMomentum(iX,iY,iZ, param2.getMomentum(iX,iY,iZ) + tau2*(-force1));
                    */
                }

                Array<T, 3> j;

                if (param1.flag(iX, iY, iZ) == interface && isWet(param2.flag(iX, iY, iZ))) {
                    // Compute the average pressure from fluid cells of both fluids.
                    T averageDensity = T();
                    plint numNeighbors = 0;
                    for (plint i = -1; i < 2; i++) {
                        for (plint j = -1; j < 2; j++) {
                            for (plint k = -1; k < 2; k++) {
                                if (isFullWet(param1.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity += param1.getDensity(iX + i, iY + j, iZ + k);
                                }
                                if (isFullWet(param2.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity +=
                                        convertRho2to1 * param2.getDensity(iX + i, iY + j, iZ + k);
                                }
                            }
                        }
                    }
                    if (numNeighbors == 0) {
                        averageDensity = rhoDefault1;
                    } else {
                        averageDensity /= (T)numNeighbors;
                    }
                    T oldDensity = param1.getDensity(iX, iY, iZ);
                    param1.setDensity(iX, iY, iZ, averageDensity);
                    param1.volumeFraction(iX, iY, iZ) = param1.mass(iX, iY, iZ) / averageDensity;
                    j = param1.getMomentum(iX, iY, iZ);
                    j *= averageDensity / oldDensity;
                    param1.setMomentum(iX, iY, iZ, j);

                    /*
                    // Implement the repelling force in terms of a momentum correction. Remember
                    that fluids can have strong
                    // density differences. In order to have the units of a momentum, this
                    correction term must be
                    // multiplied by the density.
                    T VFdiff = param1.volumeFraction(iX,iY,iZ) + param2.volumeFraction(iX,iY,iZ)
                    - 1.0; Array<T,3> normalToInterface = param1.getNormal(iX, iY, iZ);
                    //Array<T,3> normalToInterface = -param2.getNormal(iX, iY, iZ);
                    if (VFdiff>(T)0.0) {
                        T dynamicOmega2 =
                    param2.cell(iX,iY,iZ).getDynamics().getDynamicParameter(dynamicParams::dynamicOmega,
                                param2.cell(iX,iY,iZ));  // In case of a Smagorinsky model.
                        T tau2 = 0.0;
                        if (!util::isZero(dynamicOmega2)) {
                            tau2 = (T) 1 / dynamicOmega2;
                        } else {
                            tau2 = (T) 1 / param2.cell(iX,iY,iZ).getDynamics().getOmega();
                        }

                        param2.setMomentum(iX,iY,iZ, param2.getMomentum(iX,iY,iZ) +
                                interactionStrength*rhoDefault2*VFdiff*tau2*normalToInterface);
                    }
                    */
                }

                if (param2.flag(iX, iY, iZ) == interface && isWet(param1.flag(iX, iY, iZ))) {
                    // Compute the average pressure from fluid cells of both fluids.
                    T averageDensity = T();
                    plint numNeighbors = 0;
                    for (plint i = -1; i < 2; i++) {
                        for (plint j = -1; j < 2; j++) {
                            for (plint k = -1; k < 2; k++) {
                                if (isFullWet(param1.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity +=
                                        convertRho1to2 * param1.getDensity(iX + i, iY + j, iZ + k);
                                }
                                if (isFullWet(param2.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity += param2.getDensity(iX + i, iY + j, iZ + k);
                                }
                            }
                        }
                    }
                    if (numNeighbors == 0) {
                        averageDensity = rhoDefault2;
                    } else {
                        averageDensity /= (T)numNeighbors;
                    }

                    T oldDensity = param2.getDensity(iX, iY, iZ);
                    param2.setDensity(iX, iY, iZ, averageDensity);
                    param2.volumeFraction(iX, iY, iZ) = param2.mass(iX, iY, iZ) / averageDensity;
                    j = param2.getMomentum(iX, iY, iZ);
                    j *= averageDensity / oldDensity;
                    param2.setMomentum(iX, iY, iZ, j);

                    /*
                    // Implement the repelling force in terms of a momentum correction. Remember
                    that fluids can have strong
                    // density differences. In order to have the units of a momentum, this
                    correction term must be
                    // multiplied by the density.
                    T VFdiff = param1.volumeFraction(iX,iY,iZ) + param2.volumeFraction(iX,iY,iZ)
                    - 1.0; Array<T,3> normalToInterface = param2.getNormal(iX,iY,iZ);
                    //Array<T,3> normalToInterface = -param1.getNormal(iX,iY,iZ);
                    if (VFdiff>(T)0.0) {
                        T dynamicOmega1 =
                    param1.cell(iX,iY,iZ).getDynamics().getDynamicParameter(dynamicParams::dynamicOmega,
                                param1.cell(iX,iY,iZ));  // In case of a Smagorinsky model.
                        T tau1 = 0.0;
                        if (!util::isZero(dynamicOmega1)) {
                            tau1 = (T) 1 / dynamicOmega1;
                        } else {
                            tau1 = (T) 1 / param1.cell(iX,iY,iZ).getDynamics().getOmega();
                        }

                        param1.setMomentum(iX,iY,iZ, param1.getMomentum(iX,iY,iZ) +
                                interactionStrength*rhoDefault1*VFdiff*tau1*normalToInterface);
                    }
                    */
                }
            }
        }
    }
}

/* *************** Class MultiFreeSurfaceComplexCoupling3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
T MultiFreeSurfaceComplexCoupling3D<T, Descriptor>::deltaFunction(T r, T h)
{
    PLB_ASSERT(util::greaterEqual(r, (T)0));
    PLB_ASSERT(util::greaterThan(h, (T)0));

    static T pi = std::acos((T)-1);

    T delta = 0.0;
    if (util::lessEqual(r, h)) {
        delta = 0.5 * (1.0 + std::cos(pi * r / h));
    }

    return delta;
}

template <typename T, template <typename U> class Descriptor>
void MultiFreeSurfaceComplexCoupling3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    PLB_ASSERT(atomicBlocks.size() == 20);

    using namespace freeSurfaceFlag;

    std::vector<AtomicBlock3D *> atomicBlocks1, atomicBlocks2;

    for (int i = 0; i < 10; i++)
        atomicBlocks1.push_back(atomicBlocks[i]);

    for (int i = 0; i < 10; i++)
        atomicBlocks2.push_back(atomicBlocks[10 + i]);

    FreeSurfaceProcessorParam3D<T, Descriptor> param1(atomicBlocks1);
    FreeSurfaceProcessorParam3D<T, Descriptor> param2(atomicBlocks2);

    T convertRho1to2 = rhoDefault2 / rhoDefault1;
    T convertRho2to1 = rhoDefault1 / rhoDefault2;

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 3> velocity((T)0.0, (T)0.0, (T)0.0);

                if (param1.flag(iX, iY, iZ) == interface && param2.flag(iX, iY, iZ) == interface) {
                    T rhoBar1;
                    Array<T, 3> j1;
                    momentTemplates<T, Descriptor>::get_rhoBar_j(
                        param1.cell(iX, iY, iZ), rhoBar1, j1);
                    T density1 = Descriptor<T>::fullRho(rhoBar1);

                    T rhoBar2;
                    Array<T, 3> j2;
                    momentTemplates<T, Descriptor>::get_rhoBar_j(
                        param2.cell(iX, iY, iZ), rhoBar2, j2);
                    T density2 = Descriptor<T>::fullRho(rhoBar2);

                    velocity = 0.5 * (j1 / density1 + j2 / density2);
                }

                // Compute the repelling force.

                Array<T, 3> force1((T)0.0, (T)0.0, (T)0.0);

                if (isWet(param1.flag(iX, iY, iZ)) && isWet(param2.flag(iX, iY, iZ))) {
                    // Remember the envelope is 3.
                    for (plint i = -2; i < 3; i++) {
                        plint I = iX + i;
                        for (plint j = -2; j < 3; j++) {
                            plint J = iY + j;
                            for (plint k = -2; k < 3; k++) {
                                plint K = iZ + k;
                                if (isWet(param1.flag(I, J, K)) && isWet(param2.flag(I, J, K))) {
                                    T VFdiff = param1.volumeFraction(I, J, K)
                                               + param2.volumeFraction(I, J, K) - 1.0;
                                    if (VFdiff > (T)0.0) {
                                        Array<T, 3> outwardNormalToFluid1 =
                                            param1.getNormal(I, J, K) - param2.getNormal(I, J, K);
                                        outwardNormalToFluid1 /= norm(outwardNormalToFluid1);
                                        T r = std::sqrt((T)i * i + (T)j * j + (T)k * k);
                                        T delta = deltaFunction(
                                            r, (T)3);  // h is 3 because envelope is 3.
                                        force1 += -VFdiff * delta * outwardNormalToFluid1;
                                    }
                                }
                            }
                        }
                    }
                    force1 *= interactionStrength;
                }

                if (param1.flag(iX, iY, iZ) == interface && isWet(param2.flag(iX, iY, iZ))) {
                    // Compute the average pressure from fluid cells of both fluids.
                    T averageDensity = T();
                    plint numNeighbors = 0;
                    for (plint i = -1; i < 2; i++) {
                        for (plint j = -1; j < 2; j++) {
                            for (plint k = -1; k < 2; k++) {
                                if (isFullWet(param1.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity += param1.getDensity(iX + i, iY + j, iZ + k);
                                }
                                if (isFullWet(param2.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity +=
                                        convertRho2to1 * param2.getDensity(iX + i, iY + j, iZ + k);
                                }
                            }
                        }
                    }
                    if (numNeighbors == 0) {
                        averageDensity = rhoDefault1;
                    } else {
                        averageDensity /= (T)numNeighbors;
                    }
                    param1.setDensity(iX, iY, iZ, averageDensity);
                    param1.volumeFraction(iX, iY, iZ) = param1.mass(iX, iY, iZ) / averageDensity;
                    if (param2.flag(iX, iY, iZ) == interface) {
                        param1.setMomentum(iX, iY, iZ, param1.getDensity(iX, iY, iZ) * velocity);
                    }
                }

                if (param2.flag(iX, iY, iZ) == interface && isWet(param1.flag(iX, iY, iZ))) {
                    // Compute the average pressure from fluid cells of both fluids.
                    T averageDensity = T();
                    plint numNeighbors = 0;
                    for (plint i = -1; i < 2; i++) {
                        for (plint j = -1; j < 2; j++) {
                            for (plint k = -1; k < 2; k++) {
                                if (isFullWet(param1.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity +=
                                        convertRho1to2 * param1.getDensity(iX + i, iY + j, iZ + k);
                                }
                                if (isFullWet(param2.flag(iX + i, iY + j, iZ + k))) {
                                    ++numNeighbors;
                                    averageDensity += param2.getDensity(iX + i, iY + j, iZ + k);
                                }
                            }
                        }
                    }
                    if (numNeighbors == 0) {
                        averageDensity = rhoDefault2;
                    } else {
                        averageDensity /= (T)numNeighbors;
                    }
                    param2.setDensity(iX, iY, iZ, averageDensity);
                    param2.volumeFraction(iX, iY, iZ) = param2.mass(iX, iY, iZ) / averageDensity;
                    if (param1.flag(iX, iY, iZ) == interface) {
                        param2.setMomentum(iX, iY, iZ, param2.getDensity(iX, iY, iZ) * velocity);
                    }
                }

                // Implement the repelling force in terms of a momentum correction. Remember that
                // fluids can have strong density differences. In order to have the units of a
                // momentum, this correction term must be multiplied by the density.

                if (isWet(param1.flag(iX, iY, iZ)) && isWet(param2.flag(iX, iY, iZ))) {
                    // Add a force to both fluids. Remember the envelope is 3.

                    T dynamicOmega1 =
                        param1.cell(iX, iY, iZ)
                            .getDynamics()
                            .getDynamicParameter(
                                dynamicParams::dynamicOmega,
                                param1.cell(iX, iY, iZ));  // In case of a Smagorinsky model.
                    T tau1 = 0.0;
                    if (!util::isZero(dynamicOmega1)) {
                        tau1 = (T)1 / dynamicOmega1;
                    } else {
                        tau1 = (T)1 / param1.cell(iX, iY, iZ).getDynamics().getOmega();
                    }

                    T dynamicOmega2 =
                        param2.cell(iX, iY, iZ)
                            .getDynamics()
                            .getDynamicParameter(
                                dynamicParams::dynamicOmega,
                                param2.cell(iX, iY, iZ));  // In case of a Smagorinsky model.
                    T tau2 = 0.0;
                    if (!util::isZero(dynamicOmega2)) {
                        tau2 = (T)1 / dynamicOmega2;
                    } else {
                        tau2 = (T)1 / param2.cell(iX, iY, iZ).getDynamics().getOmega();
                    }

                    // force1 is an acceleration.
                    param1.setMomentum(
                        iX, iY, iZ, param1.getMomentum(iX, iY, iZ) + rhoDefault1 * tau1 * force1);
                    param2.setMomentum(
                        iX, iY, iZ,
                        param2.getMomentum(iX, iY, iZ) + rhoDefault2 * tau2 * (-force1));

                    /*
                    // force1 is a force.
                    param1.setMomentum(iX,iY,iZ, param1.getMomentum(iX,iY,iZ) + tau1*force1);
                    param2.setMomentum(iX,iY,iZ, param2.getMomentum(iX,iY,iZ) + tau2*(-force1));
                    */
                }
            }
        }
    }
}

}  // namespace plb

#endif  // MULTI_FREE_SURFACE_MODEL_3D_HH
