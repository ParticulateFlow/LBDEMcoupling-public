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
 * Division of a plane for interpolation in 3D grid refinement -- header file.
 */

#include "multiGrid/domainDivision3D.h"

namespace plb {

SurfaceDivision3D::SurfaceDivision3D(Box3D const &domain_, plint boundaryWidth, plint direction_) :
    domain(domain_), bw(boundaryWidth), direction(direction_)
{ }

Box3D SurfaceDivision3D::bulk() const
{
    // plane XY
    if (direction == 2) {
        return Box3D(
            domain.x0 + bw, domain.x1 - bw, domain.y0 + bw, domain.y1 - bw, domain.z0, domain.z1);
    }
    // plane YZ
    if (direction == 0) {
        return Box3D(
            domain.x0, domain.x1, domain.y0 + bw, domain.y1 - bw, domain.z0 + bw, domain.z1 - bw);
    }
    // plane XZ
    return Box3D(
        domain.x0 + bw, domain.x1 - bw, domain.y0, domain.y1, domain.z0 + bw, domain.z1 - bw);
}

Box3D SurfaceDivision3D::edge0N() const
{
    // plane XY
    if (direction == 2) {
        return Box3D(
            domain.x0, domain.x0 + bw - 1, domain.y0 + bw, domain.y1 - bw, domain.z0, domain.z1);
    }
    // plane YZ
    if (direction == 0) {
        return Box3D();
    }
    // plane XZ
    return Box3D();
}

Box3D SurfaceDivision3D::edge0P() const
{
    // plane XY
    if (direction == 2) {
        return Box3D(
            domain.x1 - bw + 1, domain.x1, domain.y0 + bw, domain.y1 - bw, domain.z0, domain.z1);
    }
    // plane YZ
    if (direction == 0) {
        return Box3D();
    }
    // plane XZ
    return Box3D();
}

Box3D SurfaceDivision3D::edge1N() const
{
    // plane XY
    if (direction == 2) {
        return Box3D(
            domain.x0 + bw, domain.x1 - bw, domain.y0, domain.y0 + bw - 1, domain.z0, domain.z1);
    }
    // plane YZ
    if (direction == 0) {
        return Box3D();
    }
    // plane XZ
    return Box3D();
}

Box3D SurfaceDivision3D::edge1P() const
{
    // plane XY
    if (direction == 2) {
        return Box3D(
            domain.x0 + bw, domain.x1 - bw, domain.y1 - bw + 1, domain.y1, domain.z0, domain.z1);
    }
    // plane YZ
    if (direction == 0) {
        return Box3D();
    }
    // plane XZ
    return Box3D();
}

Box3D SurfaceDivision3D::cornerNN() const
{
    // plane XY
    if (direction == 2) {
        return Box3D(
            domain.x0, domain.x0 + bw - 1, domain.y0, domain.y0 + bw - 1, domain.z0, domain.z1);
    }
    // plane YZ
    if (direction == 0) {
        return Box3D();
    }
    // plane XZ
    return Box3D();
}

Box3D SurfaceDivision3D::cornerPN() const
{
    // plane XY
    if (direction == 2) {
        return Box3D(
            domain.x1 - bw + 1, domain.x1, domain.y0, domain.y0 + bw - 1, domain.z0, domain.z1);
    }
    // plane YZ
    if (direction == 1) {
        return Box3D();
    }
    // plane XZ
    return Box3D();
}

Box3D SurfaceDivision3D::cornerNP() const
{
    // plane XY
    if (direction == 2) {
        return Box3D(
            domain.x0, domain.x0 + bw - 1, domain.y1 - bw + 1, domain.y1, domain.z0, domain.z1);
    }
    // plane YZ
    if (direction == 0) {
        return Box3D();
    }
    // plane XZ
    return Box3D();
}

Box3D SurfaceDivision3D::cornerPP() const
{
    // plane XY
    if (direction == 2) {
        return Box3D(
            domain.x1 - bw + 1, domain.x1, domain.y1 - bw + 1, domain.y1, domain.z0, domain.z1);
    }
    // plane YZ
    if (direction == 0) {
        return Box3D();
    }
    // plane XZ
    return Box3D();
}

}  // namespace plb
