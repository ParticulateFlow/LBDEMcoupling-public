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
 * Sponge (absorbing) zones, to be mainly used in addition to outflow boundary conditions -- generic
 * implementation.
 */

#ifndef SPONGE_ZONES_3D_HH
#define SPONGE_ZONES_3D_HH

#include <cmath>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "boundaryCondition/spongeZones3D.h"

namespace plb {

/* ******************** class ViscositySpongeZone3D *********************************** */

template <typename T, template <typename U> class Descriptor>
ViscositySpongeZone3D<T, Descriptor>::ViscositySpongeZone3D(
    plint nx_, plint ny_, plint nz_, T bulkOmega_, Array<plint, 6> const &numSpongeCells_,
    Array<T, 6> const &translationParameters_, Array<T, 6> const &scaleParameters_) :
    nx(nx_),
    ny(ny_),
    nz(nz_),
    bulkOmega(bulkOmega_),
    numSpongeCells(numSpongeCells_),
    translationParameters(translationParameters_),
    scaleParameters(scaleParameters_),
    useTanhSpongeFunction(true)
{
    for (int i = 0; i < 6; i++) {
        if (numSpongeCells[i] > 0) {
            PLB_ASSERT(translationParameters[i] > (T)0 && translationParameters[i] < (T)1);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ViscositySpongeZone3D<T, Descriptor>::ViscositySpongeZone3D(
    plint nx_, plint ny_, plint nz_, T bulkOmega_, Array<plint, 6> const &numSpongeCells_) :
    nx(nx_),
    ny(ny_),
    nz(nz_),
    bulkOmega(bulkOmega_),
    numSpongeCells(numSpongeCells_),
    useTanhSpongeFunction(false)
{ }

template <typename T, template <typename U> class Descriptor>
void ViscositySpongeZone3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 1);
    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(lattice);

    Dot3D offset = lattice->getLocation();

    T targetOmega = (T)1;  // Here we hardcode the value of the relaxation parameter to 1.0 for
                           // stability purposes.

    T alpha = targetOmega / bulkOmega;

    plint spongePositions[6];  // Sponge positions.
    spongePositions[0] = numSpongeCells[0];
    spongePositions[1] = nx - 1 - numSpongeCells[1];
    spongePositions[2] = numSpongeCells[2];
    spongePositions[3] = ny - 1 - numSpongeCells[3];
    spongePositions[4] = numSpongeCells[4];
    spongePositions[5] = nz - 1 - numSpongeCells[5];

    if (useTanhSpongeFunction) {
        plint spongeCenters[6];  // Sponge centers.
        spongeCenters[0] = util::roundToInt(((T)1 - translationParameters[0]) * numSpongeCells[0]);
        spongeCenters[1] =
            util::roundToInt(spongePositions[1] + translationParameters[1] * numSpongeCells[1]);
        spongeCenters[2] = util::roundToInt(((T)1 - translationParameters[2]) * numSpongeCells[2]);
        spongeCenters[3] =
            util::roundToInt(spongePositions[3] + translationParameters[3] * numSpongeCells[3]);
        spongeCenters[4] = util::roundToInt(((T)1 - translationParameters[4]) * numSpongeCells[4]);
        spongeCenters[5] =
            util::roundToInt(spongePositions[5] + translationParameters[5] * numSpongeCells[5]);

        T invSigma[6];  // Inverse sponge parameters.
        for (int i = 0; i < 6; i++) {
            invSigma[i] = (T)1 / (scaleParameters[i] * (T)numSpongeCells[i]);
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint k = iZ + offset.z;

                    T spongeFunction = (T)1;
                    bool modifyOmega = false;
                    if (numSpongeCells[0] > 0 && i <= spongePositions[0]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(i - spongeCenters[0]) * invSigma[0]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[1] > 0 && i >= spongePositions[1]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(i - spongeCenters[1]) * invSigma[1]));
                        modifyOmega = true;
                    }

                    if (numSpongeCells[2] > 0 && j <= spongePositions[2]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(j - spongeCenters[2]) * invSigma[2]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[3] > 0 && j >= spongePositions[3]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(j - spongeCenters[3]) * invSigma[3]));
                        modifyOmega = true;
                    }

                    if (numSpongeCells[4] > 0 && k <= spongePositions[4]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(k - spongeCenters[4]) * invSigma[4]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[5] > 0 && k >= spongePositions[5]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(k - spongeCenters[5]) * invSigma[5]));
                        modifyOmega = true;
                    }

                    if (modifyOmega) {
                        spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                        T localOmega = bulkOmega * spongeFunction;
                        lattice->get(iX, iY, iZ).getDynamics().setOmega(localOmega);
                    }
                }
            }
        }
    } else {
        T pi = std::acos((T)-1);
        T invNumSpongeCells[6];
        for (int i = 0; i < 6; i++) {
            invNumSpongeCells[i] = (T)1 / (T)numSpongeCells[i];
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint k = iZ + offset.z;

                    T spongeFunction = (T)1;
                    bool modifyOmega = false;
                    if (numSpongeCells[0] > 0 && i <= spongePositions[0]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)i * invNumSpongeCells[0]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[1] > 0 && i >= spongePositions[1]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(i - spongePositions[1]) * invNumSpongeCells[1]));
                        modifyOmega = true;
                    }

                    if (numSpongeCells[2] > 0 && j <= spongePositions[2]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)j * invNumSpongeCells[2]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[3] > 0 && j >= spongePositions[3]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(j - spongePositions[3]) * invNumSpongeCells[3]));
                        modifyOmega = true;
                    }

                    if (numSpongeCells[4] > 0 && k <= spongePositions[4]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)k * invNumSpongeCells[4]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[5] > 0 && k >= spongePositions[5]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(k - spongePositions[5]) * invNumSpongeCells[5]));
                        modifyOmega = true;
                    }

                    if (modifyOmega) {
                        spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                        T localOmega = bulkOmega * spongeFunction;
                        lattice->get(iX, iY, iZ).getDynamics().setOmega(localOmega);
                    }
                }
            }
        }
    }
}

/* ******************** class LocalViscositySpongeZone3D *********************************** */

template <typename T, template <typename U> class Descriptor>
LocalViscositySpongeZone3D<T, Descriptor>::LocalViscositySpongeZone3D(
    Box3D globalDomain_, T bulkOmega_, Array<plint, 6> const &numSpongeCells_) :
    globalDomain(globalDomain_), bulkOmega(bulkOmega_), numSpongeCells(numSpongeCells_)
{ }

template <typename T, template <typename U> class Descriptor>
void LocalViscositySpongeZone3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 1);
    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(lattice);

    Dot3D offset = lattice->getLocation();

    T targetOmega = (T)1;  // Here we hardcode the value of the relaxation parameter to 1.0 for
                           // stability purposes.

    T alpha = targetOmega / bulkOmega;

    plint x0 = globalDomain.x0;
    plint x1 = globalDomain.x1;
    plint y0 = globalDomain.y0;
    plint y1 = globalDomain.y1;
    plint z0 = globalDomain.z0;
    plint z1 = globalDomain.z1;

    plint spongePositions[6];  // Sponge positions.
    spongePositions[0] = x0 + numSpongeCells[0];
    spongePositions[1] = x1 - numSpongeCells[1];
    spongePositions[2] = y0 + numSpongeCells[2];
    spongePositions[3] = y1 - numSpongeCells[3];
    spongePositions[4] = z0 + numSpongeCells[4];
    spongePositions[5] = z1 - numSpongeCells[5];

    T pi = std::acos((T)-1);
    T invNumSpongeCells[6];
    for (int i = 0; i < 6; i++) {
        invNumSpongeCells[i] = (T)1 / (T)numSpongeCells[i];
    }

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint i = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint j = iY + offset.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint k = iZ + offset.z;

                T spongeFunction = (T)1;
                bool modifyOmega = false;
                if (numSpongeCells[0] > 0 && i >= x0 && i <= spongePositions[0]) {
                    spongeFunction *=
                        (T)0.5 * ((T)1 - std::cos(pi * (T)(i - x0) * invNumSpongeCells[0]));
                    modifyOmega = true;
                }
                if (numSpongeCells[1] > 0 && i >= spongePositions[1] && i <= x1) {
                    spongeFunction *=
                        (T)0.5
                        * ((T)1
                           + std::cos(pi * (T)(i - spongePositions[1]) * invNumSpongeCells[1]));
                    modifyOmega = true;
                }

                if (numSpongeCells[2] > 0 && j >= y0 && j <= spongePositions[2]) {
                    spongeFunction *=
                        (T)0.5 * ((T)1 - std::cos(pi * (T)(j - y0) * invNumSpongeCells[2]));
                    modifyOmega = true;
                }
                if (numSpongeCells[3] > 0 && j >= spongePositions[3] && j <= y1) {
                    spongeFunction *=
                        (T)0.5
                        * ((T)1
                           + std::cos(pi * (T)(j - spongePositions[3]) * invNumSpongeCells[3]));
                    modifyOmega = true;
                }

                if (numSpongeCells[4] > 0 && k >= z0 && k <= spongePositions[4]) {
                    spongeFunction *=
                        (T)0.5 * ((T)1 - std::cos(pi * (T)(k - z0) * invNumSpongeCells[4]));
                    modifyOmega = true;
                }
                if (numSpongeCells[5] > 0 && k >= spongePositions[5] && k <= z1) {
                    spongeFunction *=
                        (T)0.5
                        * ((T)1
                           + std::cos(pi * (T)(k - spongePositions[5]) * invNumSpongeCells[5]));
                    modifyOmega = true;
                }

                if (modifyOmega) {
                    spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                    T localOmega = bulkOmega * spongeFunction;
                    lattice->get(iX, iY, iZ).getDynamics().setOmega(localOmega);
                }
            }
        }
    }
}

/* ******************** class MaskedViscositySpongeZone3D *********************************** */

template <typename T, template <typename U> class Descriptor>
MaskedViscositySpongeZone3D<T, Descriptor>::MaskedViscositySpongeZone3D(
    plint nx_, plint ny_, plint nz_, T bulkOmega_, int flag_,
    Array<plint, 6> const &numSpongeCells_, Array<T, 6> const &translationParameters_,
    Array<T, 6> const &scaleParameters_) :
    nx(nx_),
    ny(ny_),
    nz(nz_),
    bulkOmega(bulkOmega_),
    numSpongeCells(numSpongeCells_),
    translationParameters(translationParameters_),
    scaleParameters(scaleParameters_),
    useTanhSpongeFunction(true),
    flag(flag_)
{
    for (int i = 0; i < 6; i++) {
        if (numSpongeCells[i] > 0) {
            PLB_ASSERT(translationParameters[i] > (T)0 && translationParameters[i] < (T)1);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
MaskedViscositySpongeZone3D<T, Descriptor>::MaskedViscositySpongeZone3D(
    plint nx_, plint ny_, plint nz_, T bulkOmega_, int flag_,
    Array<plint, 6> const &numSpongeCells_) :
    nx(nx_),
    ny(ny_),
    nz(nz_),
    bulkOmega(bulkOmega_),
    numSpongeCells(numSpongeCells_),
    useTanhSpongeFunction(false),
    flag(flag_)
{ }

template <typename T, template <typename U> class Descriptor>
void MaskedViscositySpongeZone3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 2);
    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    ScalarField3D<int> *flagMatrix = dynamic_cast<ScalarField3D<int> *>(blocks[1]);
    PLB_ASSERT(lattice);
    PLB_ASSERT(flagMatrix);

    Dot3D ofsFM = computeRelativeDisplacement(*lattice, *flagMatrix);
    Dot3D offset = lattice->getLocation();

    T targetOmega = (T)1;  // Here we hardcode the value of the relaxation parameter to 1.0 for
                           // stability purposes.

    T alpha = targetOmega / bulkOmega;

    plint spongePositions[6];  // Sponge positions.
    spongePositions[0] = numSpongeCells[0];
    spongePositions[1] = nx - 1 - numSpongeCells[1];
    spongePositions[2] = numSpongeCells[2];
    spongePositions[3] = ny - 1 - numSpongeCells[3];
    spongePositions[4] = numSpongeCells[4];
    spongePositions[5] = nz - 1 - numSpongeCells[5];

    if (useTanhSpongeFunction) {
        plint spongeCenters[6];  // Sponge centers.
        spongeCenters[0] = util::roundToInt(((T)1 - translationParameters[0]) * numSpongeCells[0]);
        spongeCenters[1] =
            util::roundToInt(spongePositions[1] + translationParameters[1] * numSpongeCells[1]);
        spongeCenters[2] = util::roundToInt(((T)1 - translationParameters[2]) * numSpongeCells[2]);
        spongeCenters[3] =
            util::roundToInt(spongePositions[3] + translationParameters[3] * numSpongeCells[3]);
        spongeCenters[4] = util::roundToInt(((T)1 - translationParameters[4]) * numSpongeCells[4]);
        spongeCenters[5] =
            util::roundToInt(spongePositions[5] + translationParameters[5] * numSpongeCells[5]);

        T invSigma[6];  // Inverse sponge parameters.
        for (int i = 0; i < 6; i++) {
            invSigma[i] = (T)1 / (scaleParameters[i] * (T)numSpongeCells[i]);
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint k = iZ + offset.z;

                    if (flagMatrix->get(iX + ofsFM.x, iY + ofsFM.y, iZ + ofsFM.z) != flag) {
                        continue;
                    }

                    T spongeFunction = (T)1;
                    bool modifyOmega = false;
                    if (numSpongeCells[0] > 0 && i <= spongePositions[0]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(i - spongeCenters[0]) * invSigma[0]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[1] > 0 && i >= spongePositions[1]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(i - spongeCenters[1]) * invSigma[1]));
                        modifyOmega = true;
                    }

                    if (numSpongeCells[2] > 0 && j <= spongePositions[2]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(j - spongeCenters[2]) * invSigma[2]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[3] > 0 && j >= spongePositions[3]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(j - spongeCenters[3]) * invSigma[3]));
                        modifyOmega = true;
                    }

                    if (numSpongeCells[4] > 0 && k <= spongePositions[4]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(k - spongeCenters[4]) * invSigma[4]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[5] > 0 && k >= spongePositions[5]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(k - spongeCenters[5]) * invSigma[5]));
                        modifyOmega = true;
                    }

                    if (modifyOmega) {
                        spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                        T localOmega = bulkOmega * spongeFunction;
                        lattice->get(iX, iY, iZ).getDynamics().setOmega(localOmega);
                    }
                }
            }
        }
    } else {
        T pi = std::acos((T)-1);
        T invNumSpongeCells[6];
        for (int i = 0; i < 6; i++) {
            invNumSpongeCells[i] = (T)1 / (T)numSpongeCells[i];
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint k = iZ + offset.z;

                    if (flagMatrix->get(iX + ofsFM.x, iY + ofsFM.y, iZ + ofsFM.z) != flag) {
                        continue;
                    }

                    T spongeFunction = (T)1;
                    bool modifyOmega = false;
                    if (numSpongeCells[0] > 0 && i <= spongePositions[0]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)i * invNumSpongeCells[0]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[1] > 0 && i >= spongePositions[1]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(i - spongePositions[1]) * invNumSpongeCells[1]));
                        modifyOmega = true;
                    }

                    if (numSpongeCells[2] > 0 && j <= spongePositions[2]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)j * invNumSpongeCells[2]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[3] > 0 && j >= spongePositions[3]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(j - spongePositions[3]) * invNumSpongeCells[3]));
                        modifyOmega = true;
                    }

                    if (numSpongeCells[4] > 0 && k <= spongePositions[4]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)k * invNumSpongeCells[4]));
                        modifyOmega = true;
                    }
                    if (numSpongeCells[5] > 0 && k >= spongePositions[5]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(k - spongePositions[5]) * invNumSpongeCells[5]));
                        modifyOmega = true;
                    }

                    if (modifyOmega) {
                        spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                        T localOmega = bulkOmega * spongeFunction;
                        lattice->get(iX, iY, iZ).getDynamics().setOmega(localOmega);
                    }
                }
            }
        }
    }
}

/* ******************** class SmagorinskySpongeZone3D *********************************** */

template <typename T, template <typename U> class Descriptor>
SmagorinskySpongeZone3D<T, Descriptor>::SmagorinskySpongeZone3D(
    plint nx_, plint ny_, plint nz_, T bulkCSmago_, T targetCSmago_,
    Array<plint, 6> const &numSpongeCells_, Array<T, 6> const &translationParameters_,
    Array<T, 6> const &scaleParameters_) :
    nx(nx_),
    ny(ny_),
    nz(nz_),
    bulkCSmago(bulkCSmago_),
    targetCSmago(targetCSmago_),
    numSpongeCells(numSpongeCells_),
    translationParameters(translationParameters_),
    scaleParameters(scaleParameters_),
    useTanhSpongeFunction(true)
{
    for (int i = 0; i < 6; i++) {
        if (numSpongeCells[i] > 0) {
            PLB_ASSERT(translationParameters[i] > (T)0 && translationParameters[i] < (T)1);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
SmagorinskySpongeZone3D<T, Descriptor>::SmagorinskySpongeZone3D(
    plint nx_, plint ny_, plint nz_, T bulkCSmago_, T targetCSmago_,
    Array<plint, 6> const &numSpongeCells_) :
    nx(nx_),
    ny(ny_),
    nz(nz_),
    bulkCSmago(bulkCSmago_),
    targetCSmago(targetCSmago_),
    numSpongeCells(numSpongeCells_),
    useTanhSpongeFunction(false)
{ }

template <typename T, template <typename U> class Descriptor>
void SmagorinskySpongeZone3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 1);
    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(lattice);

    Dot3D offset = lattice->getLocation();
    plint whichParameter = dynamicParams::smagorinskyConstant;

    T alpha = targetCSmago / bulkCSmago;

    plint spongePositions[6];  // Sponge positions.
    spongePositions[0] = numSpongeCells[0];
    spongePositions[1] = nx - 1 - numSpongeCells[1];
    spongePositions[2] = numSpongeCells[2];
    spongePositions[3] = ny - 1 - numSpongeCells[3];
    spongePositions[4] = numSpongeCells[4];
    spongePositions[5] = nz - 1 - numSpongeCells[5];

    if (useTanhSpongeFunction) {
        plint spongeCenters[6];  // Sponge centers.
        spongeCenters[0] = util::roundToInt(((T)1 - translationParameters[0]) * numSpongeCells[0]);
        spongeCenters[1] =
            util::roundToInt(spongePositions[1] + translationParameters[1] * numSpongeCells[1]);
        spongeCenters[2] = util::roundToInt(((T)1 - translationParameters[2]) * numSpongeCells[2]);
        spongeCenters[3] =
            util::roundToInt(spongePositions[3] + translationParameters[3] * numSpongeCells[3]);
        spongeCenters[4] = util::roundToInt(((T)1 - translationParameters[4]) * numSpongeCells[4]);
        spongeCenters[5] =
            util::roundToInt(spongePositions[5] + translationParameters[5] * numSpongeCells[5]);

        T invSigma[6];  // Inverse sponge parameters.
        for (int i = 0; i < 6; i++) {
            invSigma[i] = (T)1 / (scaleParameters[i] * (T)numSpongeCells[i]);
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint k = iZ + offset.z;

                    T spongeFunction = (T)1;
                    bool modifyCSmago = false;
                    if (numSpongeCells[0] > 0 && i <= spongePositions[0]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(i - spongeCenters[0]) * invSigma[0]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[1] > 0 && i >= spongePositions[1]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(i - spongeCenters[1]) * invSigma[1]));
                        modifyCSmago = true;
                    }

                    if (numSpongeCells[2] > 0 && j <= spongePositions[2]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(j - spongeCenters[2]) * invSigma[2]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[3] > 0 && j >= spongePositions[3]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(j - spongeCenters[3]) * invSigma[3]));
                        modifyCSmago = true;
                    }

                    if (numSpongeCells[4] > 0 && k <= spongePositions[4]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(k - spongeCenters[4]) * invSigma[4]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[5] > 0 && k >= spongePositions[5]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(k - spongeCenters[5]) * invSigma[5]));
                        modifyCSmago = true;
                    }

                    if (modifyCSmago) {
                        spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                        T localCSmago = bulkCSmago * spongeFunction;
                        lattice->get(iX, iY, iZ)
                            .getDynamics()
                            .setParameter(whichParameter, localCSmago);
                    }
                }
            }
        }
    } else {
        T pi = std::acos((T)-1);
        T invNumSpongeCells[6];
        for (int i = 0; i < 6; i++) {
            invNumSpongeCells[i] = (T)1 / (T)numSpongeCells[i];
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint k = iZ + offset.z;

                    T spongeFunction = (T)1;
                    bool modifyCSmago = false;
                    if (numSpongeCells[0] > 0 && i <= spongePositions[0]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)i * invNumSpongeCells[0]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[1] > 0 && i >= spongePositions[1]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(i - spongePositions[1]) * invNumSpongeCells[1]));
                        modifyCSmago = true;
                    }

                    if (numSpongeCells[2] > 0 && j <= spongePositions[2]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)j * invNumSpongeCells[2]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[3] > 0 && j >= spongePositions[3]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(j - spongePositions[3]) * invNumSpongeCells[3]));
                        modifyCSmago = true;
                    }

                    if (numSpongeCells[4] > 0 && k <= spongePositions[4]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)k * invNumSpongeCells[4]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[5] > 0 && k >= spongePositions[5]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(k - spongePositions[5]) * invNumSpongeCells[5]));
                        modifyCSmago = true;
                    }

                    // Scale and translate the sponge function.
                    if (modifyCSmago) {
                        spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                        T localCSmago = bulkCSmago * spongeFunction;
                        lattice->get(iX, iY, iZ)
                            .getDynamics()
                            .setParameter(whichParameter, localCSmago);
                    }
                }
            }
        }
    }
}

/* ******************** class MaskedSmagorinskySpongeZone3D *********************************** */

template <typename T, template <typename U> class Descriptor>
MaskedSmagorinskySpongeZone3D<T, Descriptor>::MaskedSmagorinskySpongeZone3D(
    plint nx_, plint ny_, plint nz_, T bulkCSmago_, T targetCSmago_, int flag_,
    Array<plint, 6> const &numSpongeCells_, Array<T, 6> const &translationParameters_,
    Array<T, 6> const &scaleParameters_) :
    nx(nx_),
    ny(ny_),
    nz(nz_),
    bulkCSmago(bulkCSmago_),
    targetCSmago(targetCSmago_),
    numSpongeCells(numSpongeCells_),
    translationParameters(translationParameters_),
    scaleParameters(scaleParameters_),
    useTanhSpongeFunction(true),
    flag(flag_)
{
    for (int i = 0; i < 6; i++) {
        if (numSpongeCells[i] > 0) {
            PLB_ASSERT(translationParameters[i] > (T)0 && translationParameters[i] < (T)1);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
MaskedSmagorinskySpongeZone3D<T, Descriptor>::MaskedSmagorinskySpongeZone3D(
    plint nx_, plint ny_, plint nz_, T bulkCSmago_, T targetCSmago_, int flag_,
    Array<plint, 6> const &numSpongeCells_) :
    nx(nx_),
    ny(ny_),
    nz(nz_),
    bulkCSmago(bulkCSmago_),
    targetCSmago(targetCSmago_),
    numSpongeCells(numSpongeCells_),
    useTanhSpongeFunction(false),
    flag(flag_)
{ }

template <typename T, template <typename U> class Descriptor>
void MaskedSmagorinskySpongeZone3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    PLB_ASSERT(blocks.size() == 2);
    BlockLattice3D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    ScalarField3D<int> *flagMatrix = dynamic_cast<ScalarField3D<int> *>(blocks[1]);
    PLB_ASSERT(lattice);
    PLB_ASSERT(flagMatrix);

    Dot3D ofsFM = computeRelativeDisplacement(*lattice, *flagMatrix);
    Dot3D offset = lattice->getLocation();
    plint whichParameter = dynamicParams::smagorinskyConstant;

    T alpha = targetCSmago / bulkCSmago;

    plint spongePositions[6];  // Sponge positions.
    spongePositions[0] = numSpongeCells[0];
    spongePositions[1] = nx - 1 - numSpongeCells[1];
    spongePositions[2] = numSpongeCells[2];
    spongePositions[3] = ny - 1 - numSpongeCells[3];
    spongePositions[4] = numSpongeCells[4];
    spongePositions[5] = nz - 1 - numSpongeCells[5];

    if (useTanhSpongeFunction) {
        plint spongeCenters[6];  // Sponge centers.
        spongeCenters[0] = util::roundToInt(((T)1 - translationParameters[0]) * numSpongeCells[0]);
        spongeCenters[1] =
            util::roundToInt(spongePositions[1] + translationParameters[1] * numSpongeCells[1]);
        spongeCenters[2] = util::roundToInt(((T)1 - translationParameters[2]) * numSpongeCells[2]);
        spongeCenters[3] =
            util::roundToInt(spongePositions[3] + translationParameters[3] * numSpongeCells[3]);
        spongeCenters[4] = util::roundToInt(((T)1 - translationParameters[4]) * numSpongeCells[4]);
        spongeCenters[5] =
            util::roundToInt(spongePositions[5] + translationParameters[5] * numSpongeCells[5]);

        T invSigma[6];  // Inverse sponge parameters.
        for (int i = 0; i < 6; i++) {
            invSigma[i] = (T)1 / (scaleParameters[i] * (T)numSpongeCells[i]);
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint k = iZ + offset.z;

                    if (flagMatrix->get(iX + ofsFM.x, iY + ofsFM.y, iZ + ofsFM.z) != flag) {
                        continue;
                    }

                    T spongeFunction = (T)1;
                    bool modifyCSmago = false;
                    if (numSpongeCells[0] > 0 && i <= spongePositions[0]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(i - spongeCenters[0]) * invSigma[0]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[1] > 0 && i >= spongePositions[1]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(i - spongeCenters[1]) * invSigma[1]));
                        modifyCSmago = true;
                    }

                    if (numSpongeCells[2] > 0 && j <= spongePositions[2]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(j - spongeCenters[2]) * invSigma[2]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[3] > 0 && j >= spongePositions[3]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(j - spongeCenters[3]) * invSigma[3]));
                        modifyCSmago = true;
                    }

                    if (numSpongeCells[4] > 0 && k <= spongePositions[4]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 + std::tanh((T)(k - spongeCenters[4]) * invSigma[4]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[5] > 0 && k >= spongePositions[5]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::tanh((T)(k - spongeCenters[5]) * invSigma[5]));
                        modifyCSmago = true;
                    }

                    if (modifyCSmago) {
                        spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                        T localCSmago = bulkCSmago * spongeFunction;
                        lattice->get(iX, iY, iZ)
                            .getDynamics()
                            .setParameter(whichParameter, localCSmago);
                    }
                }
            }
        }
    } else {
        T pi = std::acos((T)-1);
        T invNumSpongeCells[6];
        for (int i = 0; i < 6; i++) {
            invNumSpongeCells[i] = (T)1 / (T)numSpongeCells[i];
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint k = iZ + offset.z;

                    if (flagMatrix->get(iX + ofsFM.x, iY + ofsFM.y, iZ + ofsFM.z) != flag) {
                        continue;
                    }

                    T spongeFunction = (T)1;
                    bool modifyCSmago = false;
                    if (numSpongeCells[0] > 0 && i <= spongePositions[0]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)i * invNumSpongeCells[0]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[1] > 0 && i >= spongePositions[1]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(i - spongePositions[1]) * invNumSpongeCells[1]));
                        modifyCSmago = true;
                    }

                    if (numSpongeCells[2] > 0 && j <= spongePositions[2]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)j * invNumSpongeCells[2]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[3] > 0 && j >= spongePositions[3]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(j - spongePositions[3]) * invNumSpongeCells[3]));
                        modifyCSmago = true;
                    }

                    if (numSpongeCells[4] > 0 && k <= spongePositions[4]) {
                        spongeFunction *=
                            (T)0.5 * ((T)1 - std::cos(pi * (T)k * invNumSpongeCells[4]));
                        modifyCSmago = true;
                    }
                    if (numSpongeCells[5] > 0 && k >= spongePositions[5]) {
                        spongeFunction *=
                            (T)0.5
                            * ((T)1
                               + std::cos(pi * (T)(k - spongePositions[5]) * invNumSpongeCells[5]));
                        modifyCSmago = true;
                    }

                    // Scale and translate the sponge function.
                    if (modifyCSmago) {
                        spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                        T localCSmago = bulkCSmago * spongeFunction;
                        lattice->get(iX, iY, iZ)
                            .getDynamics()
                            .setParameter(whichParameter, localCSmago);
                    }
                }
            }
        }
    }
}

}  // namespace plb

#endif  // SPONGE_ZONES_3D_HH
