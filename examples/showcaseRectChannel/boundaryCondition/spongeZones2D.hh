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

#ifndef SPONGE_ZONES_2D_HH
#define SPONGE_ZONES_2D_HH

#include <cmath>

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "boundaryCondition/spongeZones2D.h"

namespace plb {

/* ******************** class ViscositySpongeZone2D *********************************** */

template <typename T, template <typename U> class Descriptor>
ViscositySpongeZone2D<T, Descriptor>::ViscositySpongeZone2D(
    plint nx_, plint ny_, T bulkOmega_, Array<plint, 4> const &numSpongeCells_,
    Array<T, 4> const &translationParameters_, Array<T, 4> const &scaleParameters_) :
    nx(nx_),
    ny(ny_),
    bulkOmega(bulkOmega_),
    numSpongeCells(numSpongeCells_),
    translationParameters(translationParameters_),
    scaleParameters(scaleParameters_),
    useTanhSpongeFunction(true)
{
    for (int i = 0; i < 4; i++) {
        if (numSpongeCells[i] > 0) {
            PLB_ASSERT(translationParameters[i] > (T)0 && translationParameters[i] < (T)1);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ViscositySpongeZone2D<T, Descriptor>::ViscositySpongeZone2D(
    plint nx_, plint ny_, T bulkOmega_, Array<plint, 4> const &numSpongeCells_) :
    nx(nx_),
    ny(ny_),
    bulkOmega(bulkOmega_),
    numSpongeCells(numSpongeCells_),
    useTanhSpongeFunction(false)
{ }

template <typename T, template <typename U> class Descriptor>
void ViscositySpongeZone2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_ASSERT(blocks.size() == 1);
    BlockLattice2D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(lattice);

    Dot2D offset = lattice->getLocation();

    T targetOmega = (T)1;  // Here we hardcode the value of the relaxation parameter to 1.0 for
                           // stability purposes.

    T alpha = targetOmega / bulkOmega;

    plint spongePositions[4];  // Sponge positions.
    spongePositions[0] = numSpongeCells[0];
    spongePositions[1] = nx - 1 - numSpongeCells[1];
    spongePositions[2] = numSpongeCells[2];
    spongePositions[3] = ny - 1 - numSpongeCells[3];

    if (useTanhSpongeFunction) {
        plint spongeCenters[4];  // Sponge centers.
        spongeCenters[0] = util::roundToInt(((T)1 - translationParameters[0]) * numSpongeCells[0]);
        spongeCenters[1] =
            util::roundToInt(spongePositions[1] + translationParameters[1] * numSpongeCells[1]);
        spongeCenters[2] = util::roundToInt(((T)1 - translationParameters[2]) * numSpongeCells[2]);
        spongeCenters[3] =
            util::roundToInt(spongePositions[3] + translationParameters[3] * numSpongeCells[3]);

        T invSigma[4];  // Inverse sponge parameters.
        for (int i = 0; i < 4; i++) {
            invSigma[i] = (T)1 / (scaleParameters[i] * (T)numSpongeCells[i]);
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;

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

                if (modifyOmega) {
                    spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                    T localOmega = bulkOmega * spongeFunction;
                    lattice->get(iX, iY).getDynamics().setOmega(localOmega);
                }
            }
        }
    } else {
        T pi = std::acos((T)-1);
        T invNumSpongeCells[4];
        for (int i = 0; i < 4; i++) {
            invNumSpongeCells[i] = (T)1 / (T)numSpongeCells[i];
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;

                T spongeFunction = (T)1;
                bool modifyOmega = false;
                if (numSpongeCells[0] > 0 && i <= spongePositions[0]) {
                    spongeFunction *= (T)0.5 * ((T)1 - std::cos(pi * (T)i * invNumSpongeCells[0]));
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
                    spongeFunction *= (T)0.5 * ((T)1 - std::cos(pi * (T)j * invNumSpongeCells[2]));
                    modifyOmega = true;
                }
                if (numSpongeCells[3] > 0 && j >= spongePositions[3]) {
                    spongeFunction *=
                        (T)0.5
                        * ((T)1
                           + std::cos(pi * (T)(j - spongePositions[3]) * invNumSpongeCells[3]));
                    modifyOmega = true;
                }

                if (modifyOmega) {
                    spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                    T localOmega = bulkOmega * spongeFunction;
                    lattice->get(iX, iY).getDynamics().setOmega(localOmega);
                }
            }
        }
    }
}

/* ******************** class SmagorinskySpongeZone2D *********************************** */

template <typename T, template <typename U> class Descriptor>
SmagorinskySpongeZone2D<T, Descriptor>::SmagorinskySpongeZone2D(
    plint nx_, plint ny_, T bulkCSmago_, T targetCSmago_, Array<plint, 4> const &numSpongeCells_,
    Array<T, 4> const &translationParameters_, Array<T, 4> const &scaleParameters_) :
    nx(nx_),
    ny(ny_),
    bulkCSmago(bulkCSmago_),
    targetCSmago(targetCSmago_),
    numSpongeCells(numSpongeCells_),
    translationParameters(translationParameters_),
    scaleParameters(scaleParameters_),
    useTanhSpongeFunction(true)
{
    for (int i = 0; i < 4; i++) {
        if (numSpongeCells[i] > 0) {
            PLB_ASSERT(translationParameters[i] > (T)0 && translationParameters[i] < (T)1);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
SmagorinskySpongeZone2D<T, Descriptor>::SmagorinskySpongeZone2D(
    plint nx_, plint ny_, T bulkCSmago_, T targetCSmago_, Array<plint, 4> const &numSpongeCells_) :
    nx(nx_),
    ny(ny_),
    bulkCSmago(bulkCSmago_),
    targetCSmago(targetCSmago_),
    numSpongeCells(numSpongeCells_),
    useTanhSpongeFunction(false)
{ }

template <typename T, template <typename U> class Descriptor>
void SmagorinskySpongeZone2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_ASSERT(blocks.size() == 1);
    BlockLattice2D<T, Descriptor> *lattice =
        dynamic_cast<BlockLattice2D<T, Descriptor> *>(blocks[0]);
    PLB_ASSERT(lattice);

    Dot2D offset = lattice->getLocation();
    plint whichParameter = dynamicParams::smagorinskyConstant;

    T alpha = targetCSmago / bulkCSmago;

    plint spongePositions[4];  // Sponge positions.
    spongePositions[0] = numSpongeCells[0];
    spongePositions[1] = nx - 1 - numSpongeCells[1];
    spongePositions[2] = numSpongeCells[2];
    spongePositions[3] = ny - 1 - numSpongeCells[3];

    if (useTanhSpongeFunction) {
        plint spongeCenters[4];  // Sponge centers.
        spongeCenters[0] = util::roundToInt(((T)1 - translationParameters[0]) * numSpongeCells[0]);
        spongeCenters[1] =
            util::roundToInt(spongePositions[1] + translationParameters[1] * numSpongeCells[1]);
        spongeCenters[2] = util::roundToInt(((T)1 - translationParameters[2]) * numSpongeCells[2]);
        spongeCenters[3] =
            util::roundToInt(spongePositions[3] + translationParameters[3] * numSpongeCells[3]);

        T invSigma[4];  // Inverse sponge parameters.
        for (int i = 0; i < 4; i++) {
            invSigma[i] = (T)1 / (scaleParameters[i] * (T)numSpongeCells[i]);
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;

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

                if (modifyCSmago) {
                    spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                    T localCSmago = bulkCSmago * spongeFunction;
                    lattice->get(iX, iY).getDynamics().setParameter(whichParameter, localCSmago);
                }
            }
        }
    } else {
        T pi = std::acos((T)-1);
        T invNumSpongeCells[4];
        for (int i = 0; i < 4; i++) {
            invNumSpongeCells[i] = (T)1 / (T)numSpongeCells[i];
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY + offset.y;

                T spongeFunction = (T)1;
                bool modifyCSmago = false;
                if (numSpongeCells[0] > 0 && i <= spongePositions[0]) {
                    spongeFunction *= (T)0.5 * ((T)1 - std::cos(pi * (T)i * invNumSpongeCells[0]));
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
                    spongeFunction *= (T)0.5 * ((T)1 - std::cos(pi * (T)j * invNumSpongeCells[2]));
                    modifyCSmago = true;
                }
                if (numSpongeCells[3] > 0 && j >= spongePositions[3]) {
                    spongeFunction *=
                        (T)0.5
                        * ((T)1
                           + std::cos(pi * (T)(j - spongePositions[3]) * invNumSpongeCells[3]));
                    modifyCSmago = true;
                }

                // Scale and translate the sponge function.
                if (modifyCSmago) {
                    spongeFunction = ((T)1 - alpha) * spongeFunction + alpha;
                    T localCSmago = bulkCSmago * spongeFunction;
                    lattice->get(iX, iY).getDynamics().setParameter(whichParameter, localCSmago);
                }
            }
        }
    }
}

}  // namespace plb

#endif  // SPONGE_ZONES_2D_HH
