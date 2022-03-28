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
 * Dynamics classes used to implement grid refinement -- generic implementation.
 */
#ifndef GRID_REFINEMENT_DYNAMICS_HH
#define GRID_REFINEMENT_DYNAMICS_HH

#include "multiGrid/gridRefinementDynamics.h"

namespace plb {

/* *************** Class FineGridBoundaryDynamics *************************** */

template <typename T, template <typename U> class Descriptor>
int FineGridBoundaryDynamics<T, Descriptor>::id =
    meta::registerGeneralDynamics<T, Descriptor, FineGridBoundaryDynamics<T, Descriptor> >(
        "Boundary_FineGrid");

template <typename T, template <typename U> class Descriptor>
FineGridBoundaryDynamics<T, Descriptor>::FineGridBoundaryDynamics(
    Dynamics<T, Descriptor> *baseDynamics_, TimeCounter const &timeCounter_, plint numTimeSteps_,
    plint orderOfDecomposition_) :
    BoundaryCompositeDynamics<T, Descriptor>(baseDynamics_),
    timeCounter(timeCounter_),
    numTimeSteps(numTimeSteps_),
    orderOfDecomposition(orderOfDecomposition_),
    decomposedValuesT0(),
    decomposedValuesT1()
{ }

template <typename T, template <typename U> class Descriptor>
FineGridBoundaryDynamics<T, Descriptor>::FineGridBoundaryDynamics(
    HierarchicUnserializer &unserializer) :
    BoundaryCompositeDynamics<T, Descriptor>(0, false),
    timeCounter(defaultTimeCounter),
    numTimeSteps(0),
    orderOfDecomposition(0),
    decomposedValuesT0(),
    decomposedValuesT1()
{
    this->unserialize(unserializer);
}

template <typename T, template <typename U> class Descriptor>
void FineGridBoundaryDynamics<T, Descriptor>::serialize(HierarchicSerializer &serializer) const
{
    BoundaryCompositeDynamics<T, Descriptor>::serialize(serializer);

    serializer.addValue(timeCounter);
    serializer.addValue(numTimeSteps);
    serializer.addValue(orderOfDecomposition);
    serializer.addValue(decomposedValuesT0.size());
    serializer.addValues(decomposedValuesT0);
    serializer.addValue(decomposedValuesT1.size());
    serializer.addValues(decomposedValuesT1);
}

template <typename T, template <typename U> class Descriptor>
void FineGridBoundaryDynamics<T, Descriptor>::unserialize(HierarchicUnserializer &unserializer)
{
    BoundaryCompositeDynamics<T, Descriptor>::unserialize(unserializer);

    unserializer.readValue(timeCounter);
    unserializer.readValue(numTimeSteps);
    unserializer.readValue(orderOfDecomposition);

    typename std::vector<T>::size_type t0Size, t1Size;
    unserializer.readValue(t0Size);
    decomposedValuesT0.resize(t0Size);
    unserializer.readValues(decomposedValuesT0);

    unserializer.readValue(t1Size);
    decomposedValuesT1.resize(t1Size);
    unserializer.readValues(decomposedValuesT1);
}

template <typename T, template <typename U> class Descriptor>
FineGridBoundaryDynamics<T, Descriptor> *FineGridBoundaryDynamics<T, Descriptor>::clone() const
{
    return new FineGridBoundaryDynamics<T, Descriptor>(*this);
}

/** A FineGridBoundaryDynamics cannot be composed: the end-user needs
 *  access to the original object, to manipulate internal data,
 *  the variables decomposedValuesT0 and decomposedValuesT1.
 */
template <typename T, template <typename U> class Descriptor>
bool FineGridBoundaryDynamics<T, Descriptor>::isComposeable() const
{
    return false;
}

template <typename T, template <typename U> class Descriptor>
int FineGridBoundaryDynamics<T, Descriptor>::getId() const
{
    return id;
}

template <typename T, template <typename U> class Descriptor>
void FineGridBoundaryDynamics<T, Descriptor>::completePopulations(Cell<T, Descriptor> &cell) const
{
    // Get current iteration time from the lattice
    size_t relativeTime = timeCounter.getTime() % numTimeSteps;

    // At time t0, simply recompose the corresponding decomposed values ...
    if (relativeTime == 0) {
        if ((plint)decomposedValuesT0.size()
            != cell.getDynamics().numDecomposedVariables(orderOfDecomposition))
        {
        }
        PLB_ASSERT(
            (plint)decomposedValuesT0.size()
            == cell.getDynamics().numDecomposedVariables(orderOfDecomposition));
        cell.getDynamics().recompose(cell, decomposedValuesT0, orderOfDecomposition);
    }
    // ... and at all other times, interpolate value between t0 and t1
    else
    {
        PLB_ASSERT(
            (plint)decomposedValuesT0.size()
            == cell.getDynamics().numDecomposedVariables(orderOfDecomposition));
        PLB_ASSERT(
            (plint)decomposedValuesT1.size()
            == cell.getDynamics().numDecomposedVariables(orderOfDecomposition));
        T timeFraction = (T)relativeTime / (T)numTimeSteps;
        std::vector<T> decomposedValues(decomposedValuesT0.size());
        for (pluint iVal = 0; iVal < decomposedValues.size(); ++iVal) {
            decomposedValues[iVal] = fd::linearInterpolate(
                decomposedValuesT0[iVal], decomposedValuesT1[iVal], timeFraction);
        }
        // recompose everything in the cell
        cell.getDynamics().recompose(cell, decomposedValues, orderOfDecomposition);
    }
}

template <typename T, template <typename U> class Descriptor>
std::vector<T> &FineGridBoundaryDynamics<T, Descriptor>::getDecomposedValues(plint whichTime)
{
    PLB_PRECONDITION(whichTime == 0 || whichTime == 1);
    return whichTime == 0 ? decomposedValuesT0 : decomposedValuesT1;
}

template <typename T, template <typename U> class Descriptor>

std::vector<T> const &FineGridBoundaryDynamics<T, Descriptor>::getDecomposedValues(
    plint whichTime) const
{
    PLB_PRECONDITION(whichTime == 0 || whichTime == 1);
    return whichTime == 0 ? decomposedValuesT0 : decomposedValuesT1;
}

}  // namespace plb

#endif  // GRID_REFINEMENT_DYNAMICS_HH
