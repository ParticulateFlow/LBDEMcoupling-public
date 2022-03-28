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

#ifndef ADVECTION_DIFFUSION_BOUNDARIES_H
#define ADVECTION_DIFFUSION_BOUNDARIES_H

#include "boundaryCondition/boundaryDynamics.h"
#include "complexDynamics/advectionDiffusionDynamics.h"
#include "core/globalDefs.h"

namespace plb {

/// Advection-diffusion dynamics on flat boundaries
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class AdvectionDiffusionBoundaryDynamics : public StoreDensityDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision_ = true);
    AdvectionDiffusionBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Clone the object, based on its dynamic type
    virtual AdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation> *clone()
        const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
};

/// Advection-diffusion dynamics on flat boundaries
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class RegularizedAdvectionDiffusionBoundaryDynamics : public StoreDensityDynamics<T, Descriptor> {
public:
    /// Constructor
    RegularizedAdvectionDiffusionBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision_ = true);
    RegularizedAdvectionDiffusionBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Clone the object, based on its dynamic type
    virtual RegularizedAdvectionDiffusionBoundaryDynamics<T, Descriptor, direction, orientation>
        *clone() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
};

/// Advection-diffusion dynamics on flat boundaries
template <typename T, template <typename U> class Descriptor, int direction, int orientation>
class CompleteRegularizedAdvectionDiffusionBoundaryDynamics :
    public StoreDensityDynamics<T, Descriptor> {
public:
    /// Constructor
    CompleteRegularizedAdvectionDiffusionBoundaryDynamics(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision_ = true);
    CompleteRegularizedAdvectionDiffusionBoundaryDynamics(HierarchicUnserializer &unserializer);

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Clone the object, based on its dynamic type
    virtual CompleteRegularizedAdvectionDiffusionBoundaryDynamics<
        T, Descriptor, direction, orientation>
        *clone() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
};

/// Advection-diffusion dynamics on 2D corners
template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal>
class AdvectionDiffusionCornerDynamics2D : public StoreDensityDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionCornerDynamics2D(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision_ = true);
    AdvectionDiffusionCornerDynamics2D(HierarchicUnserializer &unserializer);

    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionCornerDynamics2D<T, Descriptor, xNormal, yNormal> *clone() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
};

/// Advection-diffusion dynamics on 3D edges
template <typename T, template <typename U> class Descriptor, int plane, int normal1, int normal2>
class AdvectionDiffusionEdgeDynamics3D : public StoreDensityDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionEdgeDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision_ = true);
    AdvectionDiffusionEdgeDynamics3D(HierarchicUnserializer &unserializer);

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Clone the object, based on its dynamic type
    virtual AdvectionDiffusionEdgeDynamics3D<T, Descriptor, plane, normal1, normal2> *clone() const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
};

/// Advection-diffusion dynamics on 3D corners
template <typename T, template <typename U> class Descriptor, int xNormal, int yNormal, int zNormal>
class AdvectionDiffusionCornerDynamics3D : public StoreDensityDynamics<T, Descriptor> {
public:
    /// Constructor
    AdvectionDiffusionCornerDynamics3D(
        Dynamics<T, Descriptor> *baseDynamics, bool automaticPrepareCollision_ = true);
    AdvectionDiffusionCornerDynamics3D(HierarchicUnserializer &unserializer);

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Clone the object on its dynamic type.
    virtual AdvectionDiffusionCornerDynamics3D<T, Descriptor, xNormal, yNormal, zNormal> *clone()
        const;

    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T, Descriptor> &cell) const;

private:
    static int id;
};

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_BOUNDARIES_H
