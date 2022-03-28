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
 * Dynamics classes used to implement grid refinement -- header file.
 */

#ifndef GRID_REFINEMENT_DYNAMICS_H
#define GRID_REFINEMENT_DYNAMICS_H

#include <vector>

#include "boundaryCondition/boundaryDynamics.h"
#include "core/blockStatistics.h"
#include "core/globalDefs.h"

namespace plb {

/// Dynamics attributed to boundary cell of fine lattice with grid refinement
/** This Dynamics
 *  - Stores populations at time t (variable t0)
 *  - Stores populations at time t+1 (variable t1)
 *  - During collision, interpolates between t0 and t1, and then
 *    executes ordinary collision
 */
template <typename T, template <typename U> class Descriptor>
class FineGridBoundaryDynamics : public BoundaryCompositeDynamics<T, Descriptor> {
public:
    /// Constructor
    /** \param referenceLattice_ From the reference lattice, the FineGridBoundaryDynamics
     *                           determines the value of the current iteration.
     *  \param numTimeSteps_ Number of iteration steps leading from time t0 to time t1
     */
    FineGridBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics_, TimeCounter const &timeCounter_,
        plint numTimeSteps_, plint orderOfDecomposition_);
    FineGridBoundaryDynamics(HierarchicUnserializer &unserializer);
    virtual FineGridBoundaryDynamics<T, Descriptor> *clone() const;
    bool isComposeable() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer &serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer &unserializer);
    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;
    std::vector<T> &getDecomposedValues(plint whichTime);
    std::vector<T> const &getDecomposedValues(plint whichTime) const;

private:
    TimeCounter defaultTimeCounter;
    TimeCounter const &timeCounter;
    plint numTimeSteps;
    plint orderOfDecomposition;
    std::vector<T> decomposedValuesT0;
    std::vector<T> decomposedValuesT1;

private:
    static int id;
};

}  // namespace plb

#endif  // GRID_REFINEMENT_DYNAMICS_H
