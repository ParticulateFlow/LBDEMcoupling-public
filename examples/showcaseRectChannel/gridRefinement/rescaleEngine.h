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
 * Coupling between grids of different refinement level -- header file.
 */
#ifndef GRID_REFINEMENT_HELPERS_H
#define GRID_REFINEMENT_HELPERS_H

#include <vector>

#include "core/cell.h"
#include "core/globalDefs.h"

namespace plb {

enum executionOrder {
    lattice_decomposeT1 = -1,
    lattice_decomposeT0 = -2,
    lattice_recomposeFine = -3,
    lattice_decomposeFine = -4,
    lattice_recomposeCoarse = -5,
    lattice_recomposeBc = -6,
    container_timeInterpolation = -1,
    decomposedT12_spatialInterpolation = -1,
    decomposedT12_copyAndFilter = -2,
    decomposedT12_copyBc = -3,
    decomposedT1_spatialInterpolation = -1,
    decomposedFine_copyAndFilter = -1
};

/// A policy for scaling the data between the cells of a coarse and a fine grid.
template <typename T, template <typename U> class Descriptor>
class Rescaler {
public:
    /// Constructor
    Rescaler();

    /// Destructor
    virtual ~Rescaler();

    /// Compute the relaxation frequencies for the fine grid
    /// provided the relaxation frenquencies are known
    virtual Array<T, Descriptor<T>::q> computeRescaledRelFreq(
        const Array<T, Descriptor<T>::q> &relFreq, T xDt) const;

    /// Decompose the values of a cell at order "order"
    virtual void rescale(
        const Dynamics<T, Descriptor> &dyn, T xDt, std::vector<T> &rawData) const = 0;

    /// Get a clone of this object.
    virtual Rescaler<T, Descriptor> *clone() const = 0;

    /// Decompose the values of a cell, and rescale them to the units with respect to xDt and
    /// relFreq vector.
    virtual void decomposeAndRescale(
        Cell<T, Descriptor> const &cell, T xDt, plint order,
        std::vector<T> &decompAndRescaled) const;
};

/// Rescale values in a convective regime, dx=dt, with a factor 2 between coarse and fine grid.
/// SRT stands for single relaxatiopn time
template <typename T, template <typename U> class Descriptor>
class ConvectiveNoForceRescaler : public Rescaler<T, Descriptor> {
public:
    /// Constructor
    ConvectiveNoForceRescaler();

    /// Decompose the values of a cell at order "order"
    virtual void rescale(const Dynamics<T, Descriptor> &dyn, T xDt, std::vector<T> &rawData) const;

    /// Get a clone of this object.
    virtual ConvectiveNoForceRescaler<T, Descriptor> *clone() const;
};

/// Rescale values in a convective regime, dx=dt, with a factor 2 between coarse and fine grid.
/// SRT stands for single relaxatiopn time
template <typename T, template <typename U> class Descriptor>
class ConvectiveNoForceCompleteRescaler : public Rescaler<T, Descriptor> {
public:
    /// Constructor
    ConvectiveNoForceCompleteRescaler();

    /// Decompose the values of a cell at order "order"
    virtual void rescale(const Dynamics<T, Descriptor> &dyn, T xDt, std::vector<T> &rawData) const;

    /// Get a clone of this object.
    virtual ConvectiveNoForceCompleteRescaler<T, Descriptor> *clone() const;
};

/// Rescale values in a convective regime, dx=dt, with a factor 2 between coarse and fine grid.
/// Here the rescaling also takes into account the smagorinsky subgridscale tensor terms.
template <typename T, template <typename U> class Descriptor>
class ConvectiveNoForceCompleteSmagorinskyRescaler : public Rescaler<T, Descriptor> {
public:
    /// Constructor
    ConvectiveNoForceCompleteSmagorinskyRescaler();

    /// Decompose the values of a cell at order "order"
    virtual void rescale(const Dynamics<T, Descriptor> &dyn, T xDt, std::vector<T> &rawData) const;

    /// Get a clone of this object.
    virtual ConvectiveNoForceCompleteSmagorinskyRescaler<T, Descriptor> *clone() const;
};

}  // namespace plb

#endif  // GRID_REFINEMENT_HELPERS_H
