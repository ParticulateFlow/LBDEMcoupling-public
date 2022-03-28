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
 * Sponge (absorbing) zones, to be mainly used in addition to outflow boundary conditions -- header
 * file.
 */

#ifndef SPONGE_ZONES_2D_H
#define SPONGE_ZONES_2D_H

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

// Data processor to implement a viscosity sponge zone:
// The relaxation parameter (kinematic viscosity) is progressively increased
// inside the zone.
// The dynamics object of every cell is changed by this data processor, so
// the user must make sure in the construction of the MultiBlockLattice2D that
// each node has its own dynamics object.
template <typename T, template <typename U> class Descriptor>
class ViscositySpongeZone2D : public BoxProcessingFunctional2D {
public:
    // Constructor for the tanh sponge function.
    //   Nice value for the translation parameters is 0.5.
    //   Nice value for the scale parameters is 0.12.
    ViscositySpongeZone2D(
        plint nx_, plint ny_, T bulkOmega_, Array<plint, 4> const &numSpongeCells_,
        Array<T, 4> const &translationParameters_, Array<T, 4> const &scaleParameters_);

    // Constructor for the cos sponge function.
    ViscositySpongeZone2D(
        plint nx_, plint ny_, T bulkOmega_, Array<plint, 4> const &numSpongeCells_);

    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> blocks);

    virtual ViscositySpongeZone2D<T, Descriptor> *clone() const
    {
        return new ViscositySpongeZone2D<T, Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::dynamicVariables;  // Block lattice.
    }

private:
    plint nx, ny;  // Lattice dimensions.
    T bulkOmega;   // Value of the relaxation parameter outside of the sponge zone.
    Array<plint, 4> numSpongeCells;     // Width of the sponge zones.
    Array<T, 4> translationParameters;  // Translation parameters of the tanh sponge functions.
    Array<T, 4> scaleParameters;        // Scaling parameters of the tanh sponge functions.
    bool useTanhSpongeFunction;         // Use a tanh sponge function, or a cos sponge function.
};

// Data processor to implement a Smagorinsky sponge zone:
// The Smagorinsky parameter is progressively increased
// inside the zone.
// The dynamics object of every cell is changed by this data processor, so
// the user must make sure in the construction of the MultiBlockLattice2D that
// each node has its own dynamics object.
template <typename T, template <typename U> class Descriptor>
class SmagorinskySpongeZone2D : public BoxProcessingFunctional2D {
public:
    // Constructor for the tanh sponge function.
    //   Nice value for the translation parameters is 0.5.
    //   Nice value for the scale parameters is 0.12.
    SmagorinskySpongeZone2D(
        plint nx_, plint ny_, T bulkCSmago_, T targetCSmago_,
        Array<plint, 4> const &numSpongeCells_, Array<T, 4> const &translationParameters_,
        Array<T, 4> const &scaleParameters_);

    // Constructor for the cos sponge function.
    SmagorinskySpongeZone2D(
        plint nx_, plint ny_, T bulkCSmago_, T targetCSmago_,
        Array<plint, 4> const &numSpongeCells_);

    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> blocks);

    virtual SmagorinskySpongeZone2D<T, Descriptor> *clone() const
    {
        return new SmagorinskySpongeZone2D<T, Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::dynamicVariables;  // Block lattice.
    }

private:
    plint nx, ny;                       // Lattice dimensions.
    T bulkCSmago, targetCSmago;         // Varying parameter: bulk and target values.
    Array<plint, 4> numSpongeCells;     // Width of the sponge zones.
    Array<T, 4> translationParameters;  // Translation parameters of the tanh sponge functions.
    Array<T, 4> scaleParameters;        // Scaling parameters of the tanh sponge functions.
    bool useTanhSpongeFunction;         // Use a tanh sponge function, or a cos sponge function.
};

}  // namespace plb

#endif  // SPONGE_ZONES_2D_H
