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
#ifndef GRID_REFINEMENT_H
#define GRID_REFINEMENT_H

#include <vector>

#include "core/cell.h"
#include "core/globalDefs.h"

namespace plb {

/// Computation of Lagrange polynomial in 2D for the interpolation
template <typename T>
T interpolateValue(std::vector<T> x, std::vector<T> y);

/// A policy for scaling the data between the cells of a coarse and a fine grid.
template <typename T, template <typename U> class Descriptor>
class RescaleEngine {
public:
    virtual ~RescaleEngine() { }
    /// Decompose the values of a coarse cell, and rescale them to the units of a fine cell.
    virtual void scaleCoarseFine(
        Cell<T, Descriptor> const &coarseCell, std::vector<T> &decomposedFineValues) const = 0;
    /// Decompose the values of a fine cell, and rescale them to the units of a coarse cell.
    virtual void scaleFineCoarse(
        Cell<T, Descriptor> const &fineCell, std::vector<T> &decomposedCoarseValues) const = 0;
    /// Recompose the values into the variables of a cell, without further rescaling.
    virtual void recompose(
        Cell<T, Descriptor> &cell, std::vector<T> const &decomposedValues) const = 0;
    /// Get the order (wrt the Chapman-Enskog expansion) at which decomposition is performed.
    virtual plint getDecompositionOrder() const = 0;
    /// Get a clone of this object.
    virtual RescaleEngine<T, Descriptor> *clone() const = 0;
};

/// Rescale values in a convective regime, dx=dt, with a factor 2 between coarse and fine grid.
template <typename T, template <typename U> class Descriptor>
class ConvectiveRescaleEngine : public RescaleEngine<T, Descriptor> {
public:
    ConvectiveRescaleEngine(plint order_);
    /// Decompose the values of a coarse cell, and rescale them to the units of a fine cell.
    virtual void scaleCoarseFine(
        Cell<T, Descriptor> const &coarseCell, std::vector<T> &decomposedFineValues) const;
    /// Decompose the values of a fine cell, and rescale them to the units of a coarse cell.
    virtual void scaleFineCoarse(
        Cell<T, Descriptor> const &fineCell, std::vector<T> &decomposedCoarseValues) const;
    /// Recompose the values into the variables of a cell, without further rescaling.
    virtual void recompose(Cell<T, Descriptor> &cell, std::vector<T> const &decomposedValues) const;
    /// Get the order (wrt the Chapman-Enskog expansion) at which decomposition is performed.
    virtual plint getDecompositionOrder() const;
    /// Get a clone of this object.
    virtual ConvectiveRescaleEngine<T, Descriptor> *clone() const;

private:
    plint order;  //< Order of the decomposition wrt to the Chapman-Enskog expansion.
    static const T toFine_xDxInv;    // 1/xDx for coarse->fine
    static const T toFine_xDt;       // xDt for coarse->fine
    static const T toCoarse_xDxInv;  // 1/xDx for fine->coarse
    static const T toCoarse_xDt;     // xDt for fine->coarse
};

/// Perform no rescaling.
template <typename T, template <typename U> class Descriptor>
class NoScalingEngine : public RescaleEngine<T, Descriptor> {
public:
    NoScalingEngine(plint order_);
    /// Decompose the values of a coarse cell, and rescale them to the units of a fine cell.
    virtual void scaleCoarseFine(
        Cell<T, Descriptor> const &coarseCell, std::vector<T> &decomposedFineValues) const;
    /// Decompose the values of a fine cell, and rescale them to the units of a coarse cell.
    virtual void scaleFineCoarse(
        Cell<T, Descriptor> const &fineCell, std::vector<T> &decomposedCoarseValues) const;
    /// Recompose the values into the variables of a cell, without further rescaling.
    virtual void recompose(Cell<T, Descriptor> &cell, std::vector<T> const &decomposedValues) const;
    /// Get the order (wrt the Chapman-Enskog expansion) at which decomposition is performed.
    virtual plint getDecompositionOrder() const;
    /// Get a clone of this object.
    virtual NoScalingEngine<T, Descriptor> *clone() const;

private:
    plint order;
};

template <typename T>
void linearInterpolation(
    std::vector<T> &pop1, std::vector<T> &pop2, std::vector<T> &decomposedValues);

template <typename T>
void cubicCenteredInterpolation(
    std::vector<T> &pop1, std::vector<T> &pop2, std::vector<T> &pop3, std::vector<T> &pop4,
    std::vector<T> &decomposedValues);

template <typename T>
void quadraticNonCenteredInterpolation(
    std::vector<T> &pop1, std::vector<T> &pop2, std::vector<T> &pop3,
    std::vector<T> &decomposedValues);

}  // namespace plb

#endif  // GRID_REFINEMENT_H
