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

/* Orestis Malaspinas contributed this code.
 */

#ifndef INAMURO_ANALYTICAL_DYNAMICS_H
#define INAMURO_ANALYTICAL_DYNAMICS_H

#include "boundaryCondition/boundaryDynamics.h"
#include "core/globalDefs.h"

namespace plb {

/**
 * Implementation of Inamuro velocity boundary condition following
 * the paper
 * "A non-slip boundary condition for lattice Boltzmann simulations",
 * Inamuro, Takaji; Yoshino, Masato; Ogino, Fumimaru, (1995).
 * This implementation works for the D2Q9 Lattice only.
 */
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class InamuroAnalyticalVelocityDynamics :
    public VelocityDirichletBoundaryDynamics<T, Descriptor, direction, orientation> {
public:
    /// Constructor
    InamuroAnalyticalVelocityDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision = true);
    /// Clone the object on its dynamic type.
    virtual InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
};

/**
 * Implementation of Inamuro pressure boundary condition following
 * the paper
 * "A non-slip boundary condition for lattice Boltzmann simulations",
 * Inamuro, Takaji; Yoshino, Masato; Ogino, Fumimaru, (1995).
 * This implementation works for the D2Q9 Lattice only.
 */
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class InamuroAnalyticalPressureDynamics :
    public DensityDirichletBoundaryDynamics<T, Descriptor, direction, orientation> {
public:
    /// Constructor
    InamuroAnalyticalPressureDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision = true);
    /// Clone the object on its dynamic type.
    virtual InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation> *clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
};

}  // namespace plb

#endif  // INAMURO_ANALYTICAL_DYNAMICS_H
