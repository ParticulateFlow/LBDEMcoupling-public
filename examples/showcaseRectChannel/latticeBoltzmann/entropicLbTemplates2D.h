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

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef ENTROPIC_LB_HELPERS_2D_H
#define ENTROPIC_LB_HELPERS_2D_H

#include <cmath>

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"

namespace plb {

// template<typename T>
// struct entropicLbTemplates<T, descriptors::D2Q9Descriptor>
//{

//    typedef descriptors::D2Q9Descriptor<T> Descriptor;
//    /// Computation of equilibrium distribution with an expansion
//    /// with respect to a small velocity u
//    static T equilibrium( plint iPop, T rho, Array<T,2> const& u)
//    {
//        T c_u = Descriptor::c[iPop][0]*u[0] + Descriptor::c[iPop][1]*u[1];
//        T c_u2 = c_u*c_u;
//        T c_u3 = c_u2*c_u;
//        T c_u4 = c_u3*c_u;
//        T c_u5 = c_u4*c_u;
//        T c_u6 = c_u5*c_u;
//        T c_u7 = c_u6*c_u;

//        T uSqr = u[0]*u[0] + u[1]*u[1];
//        T uSqr2 = uSqr*uSqr;
//        T uSqr3 = uSqr2*uSqr;

//        T powUx = u[0]*u[0]*u[0]*u[0]*u[0]; // u_x^5
//        T powUy = u[1]*u[1]*u[1]*u[1]*u[1]; // u_y^5

//        T C = Descriptor::c[iPop][0] * powUx + Descriptor::c[iPop][1] * powUy;

//        powUx *= u[0]; // u_x^6
//        powUy *= u[1]; // u_y^6

//        T E = powUx + powUy;

//        powUx *= u[0]; // u_x^7
//        powUy *= u[1]; // u_y^7

//        T F = Descriptor::c[iPop][0] * powUx + Descriptor::c[iPop][1] * powUy;

//        return Descriptor::t[iPop] * rho *
//              ((T)1
//                + c_u*(C*(T)81/(T)20 + uSqr2*(T)27/(T)8 - uSqr*(T)9/(T)2
//                - E*(T)81/(T)24 - uSqr3*(T)81/(T)48 + (T)3)

//                + c_u2*(uSqr2*(T)81/(T)16 - uSqr*(T)27/(T)4
//                + C*(T)243/(T)40 + (T)9/(T)2)

//                + c_u3*(uSqr2*(T)243/(T)48 - uSqr*(T)81/(T)12 + (T)27/(T)6)

//                - c_u4*uSqr*(T)243/(T)48
//                + c_u4*(T)81/(T)24

//                - c_u5*uSqr*(T)729/(T)240
//                + c_u5*(T)243/(T)120

//                + c_u6*(T)729/(T)720

//                + c_u7*(T)2187/(T)5040

//                - C*uSqr*(T)81/(T)40

//                + C*(T)27/(T)20 - uSqr3*(T)27/(T)48 - E*(T)27/(T)24
//                - F*(T)81/(T)56 - uSqr*(T)3/(T)2 + uSqr2*(T)9/(T)8
//                )
//                - Descriptor::SkordosFactor() * Descriptor::t[iPop];
//    }

//};

}  // namespace plb

#endif
