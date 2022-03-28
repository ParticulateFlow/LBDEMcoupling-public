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
#ifndef ENTROPIC_LB_HELPERS_3D_H
#define ENTROPIC_LB_HELPERS_3D_H

#include <cmath>

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"

namespace plb {

/// Whatever is commented out corresponds to an old version, but as of now
/// nobody is able to say what it is doing.

// template<typename T>
// struct entropicLbTemplates<T, descriptors::D3Q19Descriptor>
//{
//     /// Computation of equilibrium distribution with an expansion
//     /// with respect to a small velocity u
//     static T equilibrium( plint iPop, T rho, Array<T,3> const& u)
//     {
//         typedef descriptors::D3Q19Descriptor<T> L;

//        T c_u = L::c[iPop][0]*u[0] + L::c[iPop][1]*u[1] + L::c[iPop][2]*u[2];;
//        T c_u2 = c_u*c_u;
//        T c_u3 = c_u2*c_u;
//        T c_u4 = c_u3*c_u;
//        T c_u5 = c_u4*c_u;
//        T c_u6 = c_u5*c_u;
//        T c_u7 = c_u6*c_u;

//        T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
//        T uSqr2 = uSqr*uSqr;
//        T uSqr3 = uSqr2*uSqr;

//        T powUx = u[0]*u[0]*u[0]*u[0]*u[0]; // u_x^5
//        T powUy = u[1]*u[1]*u[1]*u[1]*u[1]; // u_y^5
//        T powUz = u[2]*u[2]*u[2]*u[2]*u[2]; // u_z^5

//        T C = L::c[iPop][0] * powUx + L::c[iPop][1] * powUy
//                + L::c[iPop][2] * powUz;

//        powUx *= u[0]; // u_x^6
//        powUy *= u[1]; // u_y^6
//        powUz *= u[2]; // u_z^6

//        T E = powUx + powUy + powUz;

//        powUx *= u[0]; // u_x^7
//        powUy *= u[1]; // u_y^7
//        powUz *= u[2]; // u_z^7

//        T F = L::c[iPop][0] * powUx + L::c[iPop][1] * powUy + L::c[iPop][2] * powUz;

//        return L::t[iPop] * rho *
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
//                - L::SkordosFactor()*L::t[iPop];
//    }

//};

// template<typename T>
// struct entropicLbTemplates<T, descriptors::D3Q15Descriptor>
//{
//     /// PRL 97, 010201 (2006), implemented as in Sailfish
//     static T equilibrium( plint iPop, T rho, Array<T,3> const& u)
//     {
//         typedef descriptors::D3Q15Descriptor<T> L;

//        T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
//        T uSqr2 = uSqr*uSqr;
//        T uSqr3 = uSqr2*uSqr;
//        T uSqr4 = uSqr2*uSqr2;

//        T u0sqr = u[0]*u[0];
//        T u1sqr = u[1]*u[1];
//        T u2sqr = u[2]*u[2];

//        T chi = rho * (1 -
//                     1.5 * uSqr +
//                     1.125 * uSqr2 +
//                     1.6875 * (-uSqr3 + 2 * (u1sqr + u2sqr) *
//                                    (uSqr * u0sqr + u1sqr * u2sqr) +
//                                    20 * u0sqr * u1sqr * u2sqr) +
//                     0.6328125 * uSqr4 +
//                     2.53125 * (pow(u[0],8) + pow(u[1],8) + pow(u[2],8)
//                                - 36 * u0sqr * u1sqr * u2sqr * uSqr
//                                - pow(u[0],4) * pow(u[1],4)
//                                - pow(u[1],4) * pow(u[2],4)
//                                - pow(u[0],4) * pow(u[2],4)));

//        Array<T,3> zeta;
//        for(plint i=0; i<L::d; ++i){
//            zeta[i] = 1 +
//                     3 * u[i] +
//                     4.5 * u[i]*u[i] +
//                     4.5 * u[i]*u[i]*u[i] +
//                     3.375 * pow(u[i],4);

//            // Cyclic permutation.
//            plint x = i;
//            plint y = (i + 1) % 3;
//            plint z = (i + 2) % 3;

//            zeta[i] += 3.375 * (pow(u[x],5) - 4 * u[x] * u[y]*u[y] * u[z]*u[z]);
//            zeta[i] += 5.0625* (pow(u[x],6) - 8 * u[x]*u[x] * u[y]*u[y] * u[z]*u[z]);
//            zeta[i] += 5.0625* (pow(u[x],7) - 10 * pow(u[x],3) * u[y]*u[y] * u[z]*u[z]
//                                + 2 * u[x] * u[y]*u[y] * u[z]*u[z] * uSqr);
//            zeta[i] += 1.8984375 * (pow(u[x],8) + 16 * u[x]*u[x] * u[y]*u[y] *
//                                    u[z]*u[z] * (u[y]*u[y] + u[z]*u[z]));
//        }

//        return L::t[iPop] * chi * pow(zeta[0],L::c[iPop][0]) *
//                pow(zeta[1], L::c[iPop][1]) * pow(zeta[2], L::c[iPop][2])
//                  - L::SkordosFactor()*L::t[iPop];
//    }
//};

/// This version is actually much slower than the generic verison in entropicLbTemplates.h...
// Implementation for D3Q19 equilibrium as in Sailfish
// template<typename T>
// struct entropicLbTemplates<T, descriptors::D3Q19Descriptor>
//{
//    static T equilibrium( plint iPop, T rho, Array<T,3> const& u)
//    {
//        typedef descriptors::D3Q19Descriptor<T> L;

//        T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
//        T uSqr2 = uSqr*uSqr;

//        T u0sqr = u[0]*u[0];
//        T u1sqr = u[1]*u[1];
//        T u2sqr = u[2]*u[2];

//        T chi = rho * (1 -
//                       1.5 * uSqr +
//                       1.125 * uSqr2 +
//                       1.6875 * (pow(u[0],6) + u0sqr*u0sqr * (u1sqr + u2sqr) +
//                                (u1sqr + u2sqr) * (u1sqr*u1sqr + u2sqr*u2sqr) +
//                                u0sqr * (u1sqr*u1sqr + 12 * u1sqr * u2sqr +
//                                         u2sqr*u2sqr)) +
//                       0.6328125 * (5 * pow(u[0],8) + 5 * pow(u[1],8) +
//                                    4 * pow(u[1],6) * u2sqr +
//                                    2 * u1sqr*u1sqr * u2sqr*u2sqr +
//                                    4 * u1sqr * pow(u[2],6) +
//                                    5 * pow(u[2],8) + 4 * pow(u[0],6) * (u1sqr + u2sqr) +
//                                    4 * u0sqr * (u1sqr + u2sqr) *
//                                        (u1sqr*u1sqr + 17 * u1sqr * u2sqr +
//                                         u2sqr*u2sqr) +
//                                    2 * u0sqr*u0sqr * (u1sqr*u1sqr + 36 * u1sqr *
//                                                       u2sqr + u2sqr*u2sqr)));

//        Array<T,3> zeta;
//        for(plint i=0; i<L::d; ++i){
//            zeta[i] = 1 +
//                     3 * u[i] +
//                     4.5 * u[i]*u[i] +
//                     4.5 * u[i]*u[i]*u[i] +
//                     3.375 * pow(u[i],4);

//            // Cyclic permutation.
//            plint x = i;
//            plint y = (i + 1) % 3;
//            plint z = (i + 2) % 3;

//            zeta[i] += 3.375 * (pow(u[x],5) - 2 * u[x] * u[y]*u[y] * u[z]*u[z]);
//            zeta[i] += 5.0625* (pow(u[x],6) - 4 * u[x]*u[x] * u[y]*u[y] * u[z]*u[z]);
//            zeta[i] += 5.0625* u[x] * (pow(u[x],6) + 4 * u[x]*u[x] * u[y]*u[y] * u[z]*u[z] -
//                                       u[y]*u[y] * u[z]*u[z] * (u[y]*u[y] + u[z]*u[z]));
//            zeta[i] += 1.8984375 * u[x]*u[x] * (pow(u[x],6) - 8 * u[y]*u[y] *
//                                    u[z]*u[z] * (u[y]*u[y] + u[z]*u[z]));
//        }

//        T fEq = L::t[iPop] * chi * pow(zeta[0],L::c[iPop][0]) *
//                pow(zeta[1], L::c[iPop][1]) * pow(zeta[2], L::c[iPop][2])
//                  - L::SkordosFactor()*L::t[iPop];
//        return fEq;
//    }
//};

// Ansumali, Karlin, Ottinger (2003) for D3Q19
template <typename T>
struct entropicLbTemplates<T, descriptors::D3Q19Descriptor> {
    static T equilibrium(plint iPop, T rhoBar, Array<T, 3> const &u)
    {
        typedef descriptors::D3Q19Descriptor<T> L;

        T prod = (T)1;
        for (int iD = 0; iD < 3; ++iD) {
            T u2 = u[iD] * u[iD];

            prod *= ((T)2 - sqrt(1.0 + u2 * 3.))
                    * pow((2. * u[iD] + sqrt(1.0 + u2 * 3.)) / (1.0 - u[iD]), L::c[iPop][iD]);
        }
        return ((rhoBar + L::SkordosFactor()) * prod - L::SkordosFactor()) * L::t[iPop];
    }
};
}  // namespace plb

#endif
