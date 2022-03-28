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
 * Templates for accessing external fields (and getting a zero value
 * if the external field does not exist).
 *  -- header file
 */
#ifndef EXTERNAL_FIELD_ACCESS_H
#define EXTERNAL_FIELD_ACCESS_H

#include "core/cell.h"
#include "core/globalDefs.h"

namespace plb {

// The two classes and the function which follow implement a template mechanism
//   to compute the force according to the following rule:
//   - If there is a force in the external fields (ExternalField::sizeOfForce=d), return it
//   - If there is no force in the external fields (ExternalField::sizeOfForce=0), return 0

/// Default implementation of ExternalForceAccess: return force from external scalar.
template <typename T, template <typename U> class Descriptor, plint numForceComponents>
struct ExternalForceAccess {
    static T getComponent(Cell<T, Descriptor> const &cell, plint iD)
    {
        PLB_PRECONDITION(Descriptor<T>::d == Descriptor<T>::ExternalField::sizeOfForce);
        PLB_PRECONDITION(iD < Descriptor<T>::d);
        return *(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt + iD));
    }
    static void setComponent(Cell<T, Descriptor> const &cell, plint iD, T component)
    {
        PLB_PRECONDITION(Descriptor<T>::d == Descriptor<T>::ExternalField::sizeOfForce);
        PLB_PRECONDITION(iD < Descriptor<T>::d);
        *(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt + iD)) = component;
    }
    static Array<T, numForceComponents> get(Cell<T, Descriptor> const &cell)
    {
        PLB_PRECONDITION(Descriptor<T>::d == Descriptor<T>::ExternalField::sizeOfForce);
        Array<T, numForceComponents> force;
        force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
        return force;
    }
    static void set(Cell<T, Descriptor> &cell, Array<T, numForceComponents> const &force)
    {
        PLB_PRECONDITION(Descriptor<T>::d == Descriptor<T>::ExternalField::sizeOfForce);
        force.to_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
    }
};

/// Specialization of ExternalForceAccess: return 0 if there is no external force.
template <typename T, template <typename U> class Descriptor>
struct ExternalForceAccess<T, Descriptor, 0> {
    static T getComponent(Cell<T, Descriptor> const &cell, plint iD)
    {
        return T();
    }

    static void setComponent(Cell<T, Descriptor> const &cell, plint iD, T component) { }
    static Array<T, 0> get(Cell<T, Descriptor> const &cell)
    {
        return Array<T, 0>();
    }
    static void set(Cell<T, Descriptor> &cell, Array<T, 0> const &force) { }
};

/// Automatic instantiation of ExternalForceAccess, depending on the Descriptor
template <typename T, template <typename U> class Descriptor>
T getExternalForceComponent(Cell<T, Descriptor> const &cell, plint iD)
{
    return ExternalForceAccess<
        T, Descriptor, Descriptor<T>::ExternalField::sizeOfForce>::getComponent(cell, iD);
}

template <typename T, template <typename U> class Descriptor>
void setExternalForceComponent(Cell<T, Descriptor> const &cell, plint iD, T component)
{
    ExternalForceAccess<T, Descriptor, Descriptor<T>::ExternalField::sizeOfForce>::setComponent(
        cell, iD, component);
}

template <typename T, template <typename U> class Descriptor>
Array<T, Descriptor<T>::ExternalField::sizeOfForce> getExternalForce(
    Cell<T, Descriptor> const &cell)
{
    return ExternalForceAccess<T, Descriptor, Descriptor<T>::ExternalField::sizeOfForce>::get(cell);
}

template <typename T, template <typename U> class Descriptor>
void setExternalForce(
    Cell<T, Descriptor> &cell, Array<T, Descriptor<T>::ExternalField::sizeOfForce> const &force)
{
    return ExternalForceAccess<T, Descriptor, Descriptor<T>::ExternalField::sizeOfForce>::set(
        cell, force);
}

template <typename T, template <typename U> class Descriptor, class ExternalField>
struct RhoBarJAccess {
    static T getRhoBar(Cell<T, Descriptor> const &cell)
    {
        Dynamics<T, Descriptor> const &dynamics = cell.getDynamics();
        return dynamics.computeRhoBar(cell);
    }
    static void getJ(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &j)
    {
        T rhoBar;
        Dynamics<T, Descriptor> const &dynamics = cell.getDynamics();
        dynamics.computeRhoBarJ(cell, rhoBar, j);
    }
    static void getRhoBarJ(
        Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j)
    {
        Dynamics<T, Descriptor> const &dynamics = cell.getDynamics();
        dynamics.computeRhoBarJ(cell, rhoBar, j);
    }
};

template <typename T, template <typename U> class Descriptor>
struct RhoBarJAccess<T, Descriptor, descriptors::RhoBarJdescriptor3D> {
    static void getRhoBar(Cell<T, Descriptor> const &cell)
    {
        return *(cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt));
    }
    static void getJ(Cell<T, Descriptor> const &cell, Array<T, 3> &j)
    {
        T *externalJ = cell.getExternal(Descriptor<T>::ExternalField::jBeginsAt);
        j[0] = externalJ[0];
        j[1] = externalJ[1];
        j[2] = externalJ[2];
    }
    static void getRhoBarJ(Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, 3> &j)
    {
        rhoBar = *(cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt));
        T const *externalJ = cell.getExternal(Descriptor<T>::ExternalField::jBeginsAt);
        j[0] = externalJ[0];
        j[1] = externalJ[1];
        j[2] = externalJ[2];
    }
};

template <typename T, template <typename U> class Descriptor>
struct RhoBarJAccess<T, Descriptor, descriptors::RhoBarJdescriptor2D> {
    static void getRhoBar(Cell<T, Descriptor> const &cell)
    {
        return *(cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt));
    }
    static void getJ(Cell<T, Descriptor> const &cell, Array<T, 2> &j)
    {
        T *externalJ = cell.getExternal(Descriptor<T>::ExternalField::jBeginsAt);
        j[0] = externalJ[0];
        j[1] = externalJ[1];
    }
    static void getRhoBarJ(Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, 2> &j)
    {
        rhoBar = *(cell.getExternal(Descriptor<T>::ExternalField::rhoBarBeginsAt));
        T *externalJ = cell.getExternal(Descriptor<T>::ExternalField::jBeginsAt);
        j[0] = externalJ[0];
        j[1] = externalJ[1];
    }
};

template <typename T, template <typename U> class Descriptor>
T getRhoBar(Cell<T, Descriptor> const &cell)
{
    return RhoBarJAccess<T, Descriptor, typename Descriptor<T>::ExternalField>::getRhoBar(cell);
}

template <typename T, template <typename U> class Descriptor>
void getJ(Cell<T, Descriptor> const &cell, Array<T, Descriptor<T>::d> &j)
{
    RhoBarJAccess<T, Descriptor, typename Descriptor<T>::ExternalField>::getJ(cell, j);
}

template <typename T, template <typename U> class Descriptor>
void getRhoBarJ(Cell<T, Descriptor> const &cell, T &rhoBar, Array<T, Descriptor<T>::d> &j)
{
    RhoBarJAccess<T, Descriptor, typename Descriptor<T>::ExternalField>::getRhoBarJ(
        cell, rhoBar, j);
}

}  // namespace plb

#endif  // EXTERNAL_FIELD_ACCESS_H
