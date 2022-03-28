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
 * Templates for finding indexes for a specific subset of the neighborhood
 *  -- header file
 */
#ifndef INDEX_TEMPLATES_H
#define INDEX_TEMPLATES_H

#include <algorithm>
#include <vector>

#include "core/array.h"
#include "core/globalDefs.h"

namespace plb {

namespace indexTemplates {

/// Compute the opposite of a given direction
template <typename Descriptor>
inline plint opposite(plint iPop)
{
    if (iPop == 0)
        return 0;
    if (iPop <= Descriptor::q / 2)
        return iPop + Descriptor::q / 2;
    return iPop - Descriptor::q / 2;
}

/// Compute the opposite of a given direction
template <typename Descriptor>
inline std::vector<plint> opposite(const std::vector<plint> &pops_)
{
    std::vector<plint> pops(pops_);
    for (pluint iPop = 0; iPop < pops_.size(); ++iPop) {
        pops[iPop] = opposite<Descriptor>(pops_[iPop]);
    }
    return pops;
}

template <typename Descriptor>
plint findVelocity(Array<int, Descriptor::d> const &v)
{
    for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
        bool fit = true;
        for (int iD = 0; iD < Descriptor::d; ++iD) {
            if (Descriptor::c[iPop][iD] != v[iD]) {
                fit = false;
                break;
            }
        }
        if (fit)
            return iPop;
    }
    return Descriptor::q;
}

/// Compute the index corresponding to a specular reflection
template <typename Descriptor, int orientation>
inline plint specularReflection(plint iPop)
{
    if (iPop == 0)
        return 0;
    Array<int, Descriptor::d> v;
    for (plint iD = 0; iD < Descriptor::d; ++iD)
        v[iD] = Descriptor::c[iPop][iD];
    v[orientation] = -v[orientation];

    return findVelocity<Descriptor>(v);
}

template <typename Descriptor, plint index, plint value>
class SubIndex {
private:
    SubIndex()
    {
        for (int iVel = 0; iVel < Descriptor::q; ++iVel) {
            if (Descriptor::c[iVel][index] == value) {
                indices.push_back(iVel);
            }
        }
    }

    std::vector<plint> indices;

    template <typename Descriptor_, plint index_, plint value_>
    friend std::vector<plint> const &subIndex();
};

template <typename Descriptor, plint index, plint value>
std::vector<plint> const &subIndex()
{
    static SubIndex<Descriptor, index, value> subIndexSingleton;
    return subIndexSingleton.indices;
}

/**
 * finds distributions incoming into the wall
 * but we want the ones outgoing from the wall,
 * therefore we have to take the opposite ones.
 */
template <typename Descriptor, int direction, int orientation>
class SubIndexOutgoing {
private:
    SubIndexOutgoing()  // finds the indexes outgoing from the walls
    {
        indices = indexTemplates::subIndex<Descriptor, direction, orientation>();

        for (pluint iPop = 0; iPop < indices.size(); ++iPop) {
            indices[iPop] = indexTemplates::opposite<Descriptor>(indices[iPop]);
        }
    }

    std::vector<plint> indices;

    template <typename Descriptor_, int direction_, int orientation_>
    friend std::vector<plint> const &subIndexOutgoing();
};

template <typename Descriptor, int direction, int orientation>
std::vector<plint> const &subIndexOutgoing()
{
    static SubIndexOutgoing<Descriptor, direction, orientation> subIndexOutgoingSingleton;
    return subIndexOutgoingSingleton.indices;
}

/// finds all the remaining indexes of a lattice given some other indexes
template <typename Descriptor>
std::vector<plint> remainingIndexes(const std::vector<plint> &indices)
{
    std::vector<plint> remaining;
    for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
        bool found = false;
        for (pluint jPop = 0; jPop < indices.size(); ++jPop) {
            if (indices[jPop] == iPop) {
                found = true;
            }
        }
        if (!found) {
            remaining.push_back(iPop);
        }
    }
    return remaining;
}

/// finds the indexes outgoing from a 2D corner
template <typename Descriptor, int xNormal, int yNormal>
class SubIndexIngoingCorner2D {
private:
    SubIndexIngoingCorner2D()
    {
        typedef Descriptor L;

        Array<int, L::d> vect(-xNormal, -yNormal);
        indices.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = -xNormal;
        vect[1] = 0;
        indices.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = 0;
        vect[1] = -yNormal;
        indices.push_back(indexTemplates::findVelocity<L>(vect));
    }

    std::vector<plint> indices;

    template <typename Descriptor_, int xNormal_, int yNormal_>
    friend std::vector<plint> const &subIndexIngoingCorner2D();
};

template <typename Descriptor, int xNormal, int yNormal>
std::vector<plint> const &subIndexIngoingCorner2D()
{
    static SubIndexIngoingCorner2D<Descriptor, xNormal, yNormal> subIndexIngoingCorner2DSingleton;
    return subIndexIngoingCorner2DSingleton.indices;
}

/// finds the indexes outgoing from a 2D corner
template <typename Descriptor, int xNormal, int yNormal>
class SubIndexOutgoingCorner2D {
private:
    SubIndexOutgoingCorner2D()
    {
        typedef Descriptor L;

        Array<int, L::d> vect(xNormal, yNormal);
        std::vector<plint> knownIndexes;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = xNormal;
        vect[1] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = 0;
        vect[1] = yNormal;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = 0;
        vect[1] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));

        indices = indexTemplates::remainingIndexes<L>(knownIndexes);
    }

    std::vector<plint> indices;

    template <typename Descriptor_, int xNormal_, int yNormal_>
    friend std::vector<plint> const &subIndexOutgoingCorner2D();
};

template <typename Descriptor, int xNormal, int yNormal>
std::vector<plint> const &subIndexOutgoingCorner2D()
{
    static SubIndexOutgoingCorner2D<Descriptor, xNormal, yNormal> subIndexOutgoingCorner2DSingleton;
    return subIndexOutgoingCorner2DSingleton.indices;
}

/// finds the indexes outgoing from a 2D corner
template <typename Descriptor, int xNormal, int yNormal>
class SubIndexOutgoingInternalCorner2D {
private:
    SubIndexOutgoingInternalCorner2D()
    {
        typedef Descriptor L;

        Array<int, L::d> vect(-xNormal, -yNormal);
        indices.push_back(indexTemplates::findVelocity<L>(vect));
    }

    std::vector<plint> indices;

    template <typename Descriptor_, int xNormal_, int yNormal_>
    friend std::vector<plint> const &subIndexOutgoingInternalCorner2D();
};

template <typename Descriptor, int xNormal, int yNormal>
std::vector<plint> const &subIndexOutgoingInternalCorner2D()
{
    static SubIndexOutgoingInternalCorner2D<Descriptor, xNormal, yNormal>
        subIndexOutgoingInternalCorner2DSingleton;
    return subIndexOutgoingInternalCorner2DSingleton.indices;
}

/// finds the indexes outgoing from a 3D corner
template <typename Descriptor, int xNormal, int yNormal, int zNormal>
class SubIndexOutgoingExternalCorner3D {
private:
    SubIndexOutgoingExternalCorner3D()
    {
        typedef Descriptor L;

        Array<int, L::d> vect(xNormal, yNormal, zNormal);
        std::vector<plint> knownIndexes;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));

        vect[0] = xNormal;
        vect[1] = 0;
        vect[2] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = 0;
        vect[1] = yNormal;
        vect[2] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = 0;
        vect[1] = 0;
        vect[2] = zNormal;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));

        vect[0] = xNormal;
        vect[1] = yNormal;
        vect[2] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = 0;
        vect[1] = yNormal;
        vect[2] = zNormal;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[0] = xNormal;
        vect[1] = 0;
        vect[2] = zNormal;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));

        vect[0] = 0;
        vect[1] = 0;
        vect[2] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));

        indices = indexTemplates::remainingIndexes<L>(knownIndexes);
    }

    std::vector<plint> indices;

    template <typename Descriptor_, int xNormal_, int yNormal_, int zNormal_>
    friend std::vector<plint> const &subIndexOutgoingExternalCorner3D();
};

template <typename Descriptor, int xNormal, int yNormal, int zNormal>
std::vector<plint> const &subIndexOutgoingExternalCorner3D()
{
    static SubIndexOutgoingExternalCorner3D<Descriptor, xNormal, yNormal, zNormal>
        subIndexOutgoingExternalCorner3DSingleton;
    return subIndexOutgoingExternalCorner3DSingleton.indices;
}

/// finds the indexes outgoing from a 3D edge
template <typename Descriptor, int plane, int normal1, int normal2>
class SubIndexOutgoingExternalEdge3D {
private:
    SubIndexOutgoingExternalEdge3D()
    {
        enum { direction1 = (plane + 1) % 3, direction2 = (plane + 2) % 3 };
        typedef Descriptor L;

        Array<int, L::d> vect(0, 0, 0);
        std::vector<plint> knownIndexes;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));

        vect[direction1] = 0;
        vect[direction2] = 0;
        vect[plane] = -1;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[direction1] = 0;
        vect[direction2] = 0;
        vect[plane] = +1;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));

        vect[direction1] = normal1;
        vect[direction2] = 0;
        vect[plane] = -1;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[direction1] = normal1;
        vect[direction2] = 0;
        vect[plane] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[direction1] = normal1;
        vect[direction2] = 0;
        vect[plane] = +1;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));

        vect[direction1] = 0;
        vect[direction2] = normal2;
        vect[plane] = -1;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[direction1] = 0;
        vect[direction2] = normal2;
        vect[plane] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[direction1] = 0;
        vect[direction2] = normal2;
        vect[plane] = +1;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));

        vect[direction1] = normal1;
        vect[direction2] = normal2;
        vect[plane] = -1;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[direction1] = normal1;
        vect[direction2] = normal2;
        vect[plane] = 0;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));
        vect[direction1] = normal1;
        vect[direction2] = normal2;
        vect[plane] = +1;
        knownIndexes.push_back(indexTemplates::findVelocity<L>(vect));

        indices = indexTemplates::remainingIndexes<L>(knownIndexes);
    }

    std::vector<plint> indices;

    template <typename Descriptor_, int plane_, int normal1_, int normal2_>
    friend std::vector<plint> const &subIndexOutgoingExternalEdge3D();
};

template <typename Descriptor, int plane, int normal1, int normal2>
std::vector<plint> const &subIndexOutgoingExternalEdge3D()
{
    static SubIndexOutgoingExternalEdge3D<Descriptor, plane, normal1, normal2>
        subIndexOutgoingExternalEdge3DSingleton;
    return subIndexOutgoingExternalEdge3DSingleton.indices;
}

/// finds the indexes outgoing from a 3D edge
template <typename Descriptor, int plane, int normal1, int normal2>
class SubIndexOutgoingInternalEdge3D {
private:
    SubIndexOutgoingInternalEdge3D()
    {
        enum { direction1 = (plane + 1) % 3, direction2 = (plane + 2) % 3 };

        typedef Descriptor L;

        Array<int, L::d> vect(0, 0, 0);
        vect[direction1] = -normal1;
        vect[direction2] = -normal2;
        vect[plane] = -1;
        indices.push_back(indexTemplates::findVelocity<L>(vect));
        vect[plane] = 0;
        indices.push_back(indexTemplates::findVelocity<L>(vect));
        vect[plane] = +1;
        indices.push_back(indexTemplates::findVelocity<L>(vect));
    }

    std::vector<plint> indices;

    template <typename Descriptor_, int plane_, int normal1_, int normal2_>
    friend std::vector<plint> const &subIndexOutgoingInternalEdge3D();
};

template <typename Descriptor, int plane, int normal1, int normal2>
std::vector<plint> const &subIndexOutgoingInternalEdge3D()
{
    static SubIndexOutgoingInternalEdge3D<Descriptor, plane, normal1, normal2>
        subIndexOutgoingInternalEdge3DSingleton;
    return subIndexOutgoingInternalEdge3DSingleton.indices;
}

/// finds the indexes outgoing from a 3D corner
template <typename Descriptor, int xNormal, int yNormal, int zNormal>
class SubIndexOutgoingInternalCorner3D {
private:
    SubIndexOutgoingInternalCorner3D()
    {
        typedef Descriptor L;

        Array<int, L::d> vect(-xNormal, -yNormal, -zNormal);
        indices.push_back(indexTemplates::findVelocity<L>(vect));
    }

    std::vector<plint> indices;

    template <typename Descriptor_, int xNormal_, int yNormal_, int zNormal_>
    friend std::vector<plint> const &subIndexOutgoingInternalCorner3D();
};

template <typename Descriptor, int xNormal, int yNormal, int zNormal>
std::vector<plint> const &subIndexOutgoingInternalCorner3D()
{
    static SubIndexOutgoingInternalCorner3D<Descriptor, xNormal, yNormal, zNormal>
        subIndexOutgoingInternalCorner3DSingleton;
    return subIndexOutgoingInternalCorner3DSingleton.indices;
}

}  // namespace indexTemplates

}  // namespace plb

#endif  // INDEX_TEMPLATES_H
