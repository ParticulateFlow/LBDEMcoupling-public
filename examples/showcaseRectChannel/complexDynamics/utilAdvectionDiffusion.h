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

/** \file A helper for initialising 3D boundaries -- header file.  */
#ifndef UTIL_ADVECTION_DIFFUSION_H
#define UTIL_ADVECTION_DIFFUSION_H

#include "core/globalDefs.h"

namespace plb {

namespace utilAdvDiff {

/// finds the indexes outgoing from a 2D-3D flat wall (unknowns)
template <typename Descriptor, int direction, int orientation>
class SubIndexOutgoing {
private:
    SubIndexOutgoing()
    {
        int normalX, normalY, normalZ;
        typedef Descriptor L;

        switch (direction) {
        case 0: {
            if (orientation == 1)
                normalX = 1;
            else
                normalX = -1;
            normalY = 0;
            normalZ = 0;
            break;
        }
        case 1: {
            if (orientation == 1)
                normalY = 1;
            else
                normalY = -1;
            normalX = 0;
            normalZ = 0;
            break;
        }
        case 2: {
            if (orientation == 1)
                normalZ = 1;
            else
                normalZ = -1;
            normalX = 0;
            normalY = 0;
            break;
        }
        }

        // I define the dimension of Normal Vec=3 to make this routine usable for 2D and 3D flat
        // walls
        std::vector<plint> NormalVec(3, 0);
        NormalVec[0] = normalX;
        NormalVec[1] = normalY;
        NormalVec[2] = normalZ;
        // add zero velocity
        // knownIndexes.push_back(0);
        // compute scalar product with boundary normal for all other velocities
        for (plint iPop = 1; iPop < L::q; ++iPop) {
            plint sum = 0;
            for (int id = 0; id < L::d; ++id) {
                sum += L::c[iPop][id] * NormalVec[id];
            }
            if (sum < 0) {
                indices.push_back(iPop);
            }
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

// For Egdges
template <typename Descriptor, int plane, int normal1, int normal2>
class SubIndexOutgoing3DonEdges {
private:
    SubIndexOutgoing3DonEdges()
    {
        int normalX, normalY, normalZ;
        typedef Descriptor L;

        switch (plane) {
        case 0: {
            normalX = 0;
            if (normal1 == 1)
                normalY = 1;
            else
                normalY = -1;
            if (normal2 == 1)
                normalZ = 1;
            else
                normalZ = -1;
            break;
        }
        case 1: {
            normalY = 0;
            if (normal1 == 1)
                normalZ = 1;
            else
                normalZ = -1;
            if (normal2 == 1)
                normalX = 1;
            else
                normalX = -1;
            break;
        }
        case 2: {
            normalZ = 0;
            if (normal1 == 1)
                normalX = 1;
            else
                normalX = -1;
            if (normal2 == 1)
                normalY = 1;
            else
                normalY = -1;
            break;
        }
        }

        // add zero velocity
        // knownIndexes.push_back(0);
        // compute scalar product with boundary normal for all other velocities
        for (plint iP = 1; iP < L::q; ++iP) {
            if ((L::c[iP][0] * normalX + L::c[iP][1] * normalY + L::c[iP][2] * normalZ) < 0) {
                indices.push_back(iP);
            }
        }
    }
    std::vector<plint> indices;

    template <typename Descriptor_, int plane_, int normal1_, int normal2_>
    friend std::vector<plint> const &subIndexOutgoing3DonEdges();
};

template <typename Descriptor, int plane, int normal1, int normal2>
std::vector<plint> const &subIndexOutgoing3DonEdges()
{
    static SubIndexOutgoing3DonEdges<Descriptor, plane, normal1, normal2>
        subIndexOutgoing3DonEdgesSingleton;
    return subIndexOutgoing3DonEdgesSingleton.indices;
}

// For 3D Corners
template <typename Descriptor, int normalX, int normalY, int normalZ>
class SubIndexOutgoing3DonCorners {
private:
    SubIndexOutgoing3DonCorners()
    {
        typedef Descriptor L;

        // add zero velocity
        // knownIndexes.push_back(0);
        // compute scalar product with boundary normal for all other velocities
        for (plint iP = 1; iP < L::q; ++iP) {
            if (L::c[iP][0] * normalX + L::c[iP][1] * normalY + L::c[iP][2] * normalZ < 0) {
                indices.push_back(iP);
            }
        }
    }

    std::vector<plint> indices;

    template <typename Descriptor_, int normalX_, int normalY_, int normalZ_>
    friend std::vector<plint> const &subIndexOutgoing3DonCorners();
};

template <typename Descriptor, int normalX, int normalY, int normalZ>
std::vector<plint> const &subIndexOutgoing3DonCorners()
{
    static SubIndexOutgoing3DonCorners<Descriptor, normalX, normalY, normalZ>
        subIndexOutgoing3DonCornersSingleton;
    return subIndexOutgoing3DonCornersSingleton.indices;
}

// For 2D Corners

template <typename Descriptor, int normalX, int normalY>
class SubIndexOutgoing2DonCorners {
private:
    SubIndexOutgoing2DonCorners()
    {
        typedef Descriptor L;

        // add zero velocity
        // knownIndexes.push_back(0);
        // compute scalar product with boundary normal for all other velocities
        for (plint iPop = 1; iPop < L::q; ++iPop) {
            if (L::c[iPop][0] * normalX + L::c[iPop][1] * normalY < 0) {
                indices.push_back(iPop);
            }
        }
    }

    std::vector<plint> indices;

    template <typename Descriptor_, int normalX_, int normalY_>
    friend std::vector<plint> const &subIndexOutgoing2DonCorners();
};

template <typename Descriptor, int normalX, int normalY>
std::vector<plint> const &subIndexOutgoing2DonCorners()
{
    static SubIndexOutgoing2DonCorners<Descriptor, normalX, normalY>
        subIndexOutgoing2DonCornersSingleton;
    return subIndexOutgoing2DonCornersSingleton.indices;
}

}  // namespace utilAdvDiff

}  // namespace plb

#endif
