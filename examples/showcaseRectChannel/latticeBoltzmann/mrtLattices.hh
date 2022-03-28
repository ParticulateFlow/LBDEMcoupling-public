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
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the Palabos code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- generic code
 */
#ifndef MRT_LATTICES_HH
#define MRT_LATTICES_HH

namespace plb {
namespace descriptors {

// MRT D2Q9 ////////////////////////////////////////////////////////////

template <typename T>
const T MRTD2Q9DescriptorBase<T>::M[BaseDescriptor::q][BaseDescriptor::q] = {
    {(T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1},
    {-(T)4, (T)2, -(T)1, (T)2, -(T)1, (T)2, -(T)1, (T)2, -(T)1},
    {(T)4, (T)1, -(T)2, (T)1, -(T)2, (T)1, -(T)2, (T)1, -(T)2},
    {T(), -(T)1, -(T)1, -(T)1, T(), (T)1, (T)1, (T)1, T()},
    {T(), -(T)1, (T)2, -(T)1, T(), (T)1, -(T)2, (T)1, T()},
    {T(), (T)1, T(), -(T)1, -(T)1, -(T)1, T(), (T)1, (T)1},
    {T(), (T)1, T(), -(T)1, (T)2, -(T)1, T(), (T)1, -(T)2},
    {T(), T(), (T)1, T(), -(T)1, T(), (T)1, T(), -(T)1},
    {T(), -(T)1, T(), (T)1, T(), -(T)1, T(), (T)1, T()}};

template <typename T>
const T MRTD2Q9DescriptorBase<T>::invM[BaseDescriptor::q][BaseDescriptor::q] = {
    {(T)1 / (T)9, -(T)1 / (T)9, (T)1 / (T)9, T(), T(), T(), T(), T(), T()},
    {(T)1 / (T)9, (T)1 / (T)18, (T)1 / (T)36, -(T)1 / (T)6, -(T)1 / (T)12, (T)1 / (T)6,
     (T)1 / (T)12, T(), -(T)1 / (T)4},
    {(T)1 / (T)9, -(T)1 / (T)36, -(T)1 / (T)18, -(T)1 / (T)6, (T)1 / (T)6, T(), T(), (T)1 / (T)4,
     T()},
    {(T)1 / (T)9, (T)1 / (T)18, (T)1 / (T)36, -(T)1 / (T)6, -(T)1 / (T)12, -(T)1 / (T)6,
     -(T)1 / (T)12, T(), (T)1 / (T)4},
    {(T)1 / (T)9, -(T)1 / (T)36, -(T)1 / (T)18, T(), T(), -(T)1 / (T)6, (T)1 / (T)6, -(T)1 / (T)4,
     T()},
    {(T)1 / (T)9, (T)1 / (T)18, (T)1 / (T)36, (T)1 / (T)6, (T)1 / (T)12, -(T)1 / (T)6,
     -(T)1 / (T)12, T(), -(T)1 / (T)4},
    {(T)1 / (T)9, -(T)1 / (T)36, -(T)1 / (T)18, (T)1 / (T)6, -(T)1 / (T)6, T(), T(), (T)1 / (T)4,
     T()},
    {(T)1 / (T)9, (T)1 / (T)18, (T)1 / (T)36, (T)1 / (T)6, (T)1 / (T)12, (T)1 / (T)6, (T)1 / (T)12,
     T(), (T)1 / (T)4},
    {(T)1 / (T)9, -(T)1 / (T)36, -(T)1 / (T)18, T(), T(), (T)1 / (T)6, -(T)1 / (T)6, -(T)1 / (T)4,
     T()}};

template <typename T>
const T MRTD2Q9DescriptorBase<T>::S[BaseDescriptor::q] = {(T)1, (T)1.63, (T)1.14, (T)1, (T)1.92,
                                                          (T)1, (T)1.92, T(),     T()};
// s7=s8 to have a shear viscosity nu
// and the bulk viscosity depends on s1.
// sq->s4 and s6
// s_epsilon -> s2

template <typename T>
const int MRTD2Q9DescriptorBase<T>::momentumIndexes[MRTD2Q9DescriptorBase<T>::jIndexes] = {3, 5};

template <typename T>
const int MRTD2Q9DescriptorBase<T>::shearViscIndexes[MRTD2Q9DescriptorBase<T>::shearIndexes] = {
    7, 8};

template <typename T>
const int MRTD2Q9DescriptorBase<T>::qViscIndexes[MRTD2Q9DescriptorBase<T>::qIndexes] = {4, 6};

template <typename T>
const char MRTD2Q9Descriptor<T>::name[] = "MRTD2Q9";

// EXTERNAL FORCE MRT D2Q9 ////////////////////////////////////////////////////////////

template <typename T>
const char ForcedMRTD2Q9Descriptor<T>::name[] = "ForcedMRTD2Q9";

// External Force, density and momentum for ShanChen MRT 2D model
template <typename T>
const char ForcedMRTShanChenD2Q9Descriptor<T>::name[] = "ForcedShanChenMRTD2Q9";

// EXTERNAL VELOCITY MRT D2Q9 ////////////////////////////////////////////////////////////
template <typename T>
const char ExternalVelocityMRTD2Q9Descriptor<T>::name[] = "ExternalVelocityMRTD2Q9";

// MRT D3Q15 ////////////////////////////////////////////////////////////

template <typename T>
const T MRTD3Q15DescriptorBase<T>::M[BaseDescriptor::q][BaseDescriptor::q] = {
    {(T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1},
    {-(T)2, -(T)1, -(T)1, -(T)1, (T)1, (T)1, (T)1, (T)1, -(T)1, -(T)1, -(T)1, (T)1, (T)1, (T)1,
     (T)1},
    {(T)16, -(T)4, -(T)4, -(T)4, (T)1, (T)1, (T)1, (T)1, -(T)4, -(T)4, -(T)4, (T)1, (T)1, (T)1,
     (T)1},
    {T(), -(T)1, T(), T(), -(T)1, -(T)1, -(T)1, -(T)1, (T)1, T(), T(), (T)1, (T)1, (T)1, (T)1},
    {T(), (T)4, T(), T(), -(T)1, -(T)1, -(T)1, -(T)1, -(T)4, T(), T(), (T)1, (T)1, (T)1, (T)1},
    {T(), T(), -(T)1, T(), -(T)1, -(T)1, (T)1, (T)1, T(), (T)1, T(), (T)1, (T)1, -(T)1, -(T)1},
    {T(), T(), (T)4, T(), -(T)1, -(T)1, (T)1, (T)1, T(), -(T)4, T(), (T)1, (T)1, -(T)1, -(T)1},
    {T(), T(), T(), -(T)1, -(T)1, (T)1, -(T)1, (T)1, T(), T(), (T)1, (T)1, -(T)1, (T)1, -(T)1},
    {T(), T(), T(), (T)4, -(T)1, (T)1, -(T)1, (T)1, T(), T(), -(T)4, (T)1, -(T)1, (T)1, -(T)1},
    {T(), (T)2, -(T)1, -(T)1, T(), T(), T(), T(), (T)2, -(T)1, -(T)1, T(), T(), T(), T()},
    {T(), T(), (T)1, -(T)1, T(), T(), T(), T(), T(), (T)1, -(T)1, T(), T(), T(), T()},
    {T(), T(), T(), T(), (T)1, (T)1, -(T)1, -(T)1, T(), T(), T(), (T)1, (T)1, -(T)1, -(T)1},
    {T(), T(), T(), T(), (T)1, -(T)1, -(T)1, (T)1, T(), T(), T(), (T)1, -(T)1, -(T)1, (T)1},
    {T(), T(), T(), T(), (T)1, -(T)1, (T)1, -(T)1, T(), T(), T(), (T)1, -(T)1, (T)1, -(T)1},
    {T(), T(), T(), T(), -(T)1, (T)1, (T)1, -(T)1, T(), T(), T(), (T)1, -(T)1, -(T)1, (T)1}};

template <typename T>
const T MRTD3Q15DescriptorBase<T>::invM[BaseDescriptor::q][BaseDescriptor::q] = {
    {(T)1 / (T)15, -(T)1 / (T)9, (T)2 / (T)45, T(), T(), T(), T(), T(), T(), T(), T(), T(), T(),
     T(), T()},
    {(T)1 / (T)15, -(T)1 / (T)18, -(T)1 / (T)90, -(T)1 / (T)10, (T)1 / (T)10, T(), T(), T(), T(),
     (T)1 / (T)6, T(), T(), T(), T(), T()},
    {(T)1 / (T)15, -(T)1 / (T)18, -(T)1 / (T)90, T(), T(), -(T)1 / (T)10, (T)1 / (T)10, T(), T(),
     -(T)1 / (T)12, (T)1 / (T)4, T(), T(), T(), T()},
    {(T)1 / (T)15, -(T)1 / (T)18, -(T)1 / (T)90, T(), T(), T(), T(), -(T)1 / (T)10, (T)1 / (T)10,
     -(T)1 / (T)12, -(T)1 / (T)4, T(), T(), T(), T()},
    {(T)1 / (T)15, (T)1 / (T)18, (T)1 / (T)360, -(T)1 / (T)10, -(T)1 / (T)40, -(T)1 / (T)10,
     -(T)1 / (T)40, -(T)1 / (T)10, -(T)1 / (T)40, T(), T(), (T)1 / (T)8, (T)1 / (T)8, (T)1 / (T)8,
     -(T)1 / (T)8},
    {(T)1 / (T)15, (T)1 / (T)18, (T)1 / (T)360, -(T)1 / (T)10, -(T)1 / (T)40, -(T)1 / (T)10,
     -(T)1 / (T)40, (T)1 / (T)10, (T)1 / (T)40, T(), T(), (T)1 / (T)8, -(T)1 / (T)8, -(T)1 / (T)8,
     (T)1 / (T)8},
    {(T)1 / (T)15, (T)1 / (T)18, (T)1 / (T)360, -(T)1 / (T)10, -(T)1 / (T)40, (T)1 / (T)10,
     (T)1 / (T)40, -(T)1 / (T)10, -(T)1 / (T)40, T(), T(), -(T)1 / (T)8, -(T)1 / (T)8, (T)1 / (T)8,
     (T)1 / (T)8},
    {(T)1 / (T)15, (T)1 / (T)18, (T)1 / (T)360, -(T)1 / (T)10, -(T)1 / (T)40, (T)1 / (T)10,
     (T)1 / (T)40, (T)1 / (T)10, (T)1 / (T)40, T(), T(), -(T)1 / (T)8, (T)1 / (T)8, -(T)1 / (T)8,
     -(T)1 / (T)8},
    {(T)1 / (T)15, -(T)1 / (T)18, -(T)1 / (T)90, (T)1 / (T)10, -(T)1 / (T)10, T(), T(), T(), T(),
     (T)1 / (T)6, T(), T(), T(), T(), T()},
    {(T)1 / (T)15, -(T)1 / (T)18, -(T)1 / (T)90, T(), T(), (T)1 / (T)10, -(T)1 / (T)10, T(), T(),
     -(T)1 / (T)12, (T)1 / (T)4, T(), T(), T(), T()},
    {(T)1 / (T)15, -(T)1 / (T)18, -(T)1 / (T)90, T(), T(), T(), T(), (T)1 / (T)10, -(T)1 / (T)10,
     -(T)1 / (T)12, -(T)1 / (T)4, T(), T(), T(), T()},
    {(T)1 / (T)15, (T)1 / (T)18, (T)1 / (T)360, (T)1 / (T)10, (T)1 / (T)40, (T)1 / (T)10,
     (T)1 / (T)40, (T)1 / (T)10, (T)1 / (T)40, T(), T(), (T)1 / (T)8, (T)1 / (T)8, (T)1 / (T)8,
     (T)1 / (T)8},
    {(T)1 / (T)15, (T)1 / (T)18, (T)1 / (T)360, (T)1 / (T)10, (T)1 / (T)40, (T)1 / (T)10,
     (T)1 / (T)40, -(T)1 / (T)10, -(T)1 / (T)40, T(), T(), (T)1 / (T)8, -(T)1 / (T)8, -(T)1 / (T)8,
     -(T)1 / (T)8},
    {(T)1 / (T)15, (T)1 / (T)18, (T)1 / (T)360, (T)1 / (T)10, (T)1 / (T)40, -(T)1 / (T)10,
     -(T)1 / (T)40, (T)1 / (T)10, (T)1 / (T)40, T(), T(), -(T)1 / (T)8, -(T)1 / (T)8, (T)1 / (T)8,
     -(T)1 / (T)8},
    {(T)1 / (T)15, (T)1 / (T)18, (T)1 / (T)360, (T)1 / (T)10, (T)1 / (T)40, -(T)1 / (T)10,
     -(T)1 / (T)40, -(T)1 / (T)10, -(T)1 / (T)40, T(), T(), -(T)1 / (T)8, (T)1 / (T)8, -(T)1 / (T)8,
     (T)1 / (T)8}};

template <typename T>
const int MRTD3Q15DescriptorBase<T>::momentumIndexes[MRTD3Q15DescriptorBase<T>::jIndexes] = {
    3, 5, 7};

template <typename T>
const int MRTD3Q15DescriptorBase<T>::shearViscIndexes[MRTD3Q15DescriptorBase<T>::shearIndexes] = {
    9, 10, 11, 12, 13};

template <typename T>
const int MRTD3Q15DescriptorBase<T>::qViscIndexes[MRTD3Q15DescriptorBase<T>::qIndexes] = {4, 6, 8};

template <typename T>
const T MRTD3Q15DescriptorBase<T>::S[BaseDescriptor::q] = {
    /*s0*/ (T)1,     // rho (conserved)
    /*s1*/ (T)1.1,   // bulk viscosity index
    /*s2*/ (T)1.1,   // e index
    /*s3*/ (T)1,     // rho*ux (conserved)
    /*s4*/ (T)1.1,   // <=> q index
    /*s5*/ (T)1,     // rho*uy (conserved)
    /*s6*/ (T)1.1,   // = s4 <=> q index
    /*s7*/ (T)1,     // rho*uz (conserved)
    /*s8*/ (T)1.1,   // = s4 <=> q index
    /*s9*/ T(),      // used to define nu
    /*s10*/ T(),     // = s9
    /*s11*/ T(),     // = s9,
    /*s12*/ T(),     // = s9,
    /*s13*/ T(),     // = s9,
    /*s14*/ T(),     // = s9,,
    /*s15*/ (T)1.1,  // = m index,
};

template <typename T>
const char MRTD3Q15Descriptor<T>::name[] = "MRTD3Q15";

// MRT D3Q19 ////////////////////////////////////////////////////////////

template <typename T>
const T MRTD3Q19DescriptorBase<T>::M[BaseDescriptor::q][BaseDescriptor::q] = {
    /*0*/ {
        (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1,
        (T)1, (T)1, (T)1, (T)1},
    /*1*/
    {(T)-30, (T)-11, (T)-11, (T)-11, (T)8, (T)8, (T)8, (T)8, (T)8, (T)8, (T)-11, (T)-11, (T)-11,
     (T)8, (T)8, (T)8, (T)8, (T)8, (T)8},
    /*2*/
    {(T)12, (T)-4, (T)-4, (T)-4, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)-4, (T)-4, (T)-4, (T)1,
     (T)1, (T)1, (T)1, (T)1, (T)1},
    /*3*/
    {T(), (T)-1, T(), T(), (T)-1, (T)-1, (T)-1, (T)-1, T(), T(), (T)1, T(), T(), (T)1, (T)1, (T)1,
     (T)1, T(), T()},
    /*4*/
    {T(), (T)4, T(), T(), (T)-1, (T)-1, (T)-1, (T)-1, T(), T(), (T)-4, T(), T(), (T)1, (T)1, (T)1,
     (T)1, T(), T()},
    /*5*/
    {T(), T(), (T)-1, T(), (T)-1, (T)1, T(), T(), (T)-1, (T)-1, T(), (T)1, T(), (T)1, (T)-1, T(),
     T(), (T)1, (T)1},
    /*6*/
    {T(), T(), (T)4, T(), (T)-1, (T)1, T(), T(), (T)-1, (T)-1, T(), (T)-4, T(), (T)1, (T)-1, T(),
     T(), (T)1, (T)1},
    /*7*/
    {T(), T(), T(), (T)-1, T(), T(), (T)-1, (T)1, (T)-1, (T)1, T(), T(), (T)1, T(), T(), (T)1,
     (T)-1, (T)1, (T)-1},
    /*8*/
    {T(), T(), T(), (T)4, T(), T(), (T)-1, (T)1, (T)-1, (T)1, T(), T(), (T)-4, T(), T(), (T)1,
     (T)-1, (T)1, (T)-1},
    /*9*/
    {T(), (T)2, (T)-1, (T)-1, (T)1, (T)1, (T)1, (T)1, (T)-2, (T)-2, (T)2, (T)-1, (T)-1, (T)1, (T)1,
     (T)1, (T)1, (T)-2, (T)-2},
    /*10*/
    {T(), (T)-4, (T)2, (T)2, (T)1, (T)1, (T)1, (T)1, (T)-2, (T)-2, (T)-4, (T)2, (T)2, (T)1, (T)1,
     (T)1, (T)1, (T)-2, (T)-2},
    /*11*/
    {T(), T(), (T)1, (T)-1, (T)1, (T)1, (T)-1, (T)-1, T(), T(), T(), (T)1, (T)-1, (T)1, (T)1, (T)-1,
     (T)-1, T(), T()},
    /*12*/
    {T(), T(), (T)-2, (T)2, (T)1, (T)1, (T)-1, (T)-1, T(), T(), T(), (T)-2, (T)2, (T)1, (T)1, (T)-1,
     (T)-1, T(), T()},
    /*13*/
    {T(), T(), T(), T(), (T)1, (T)-1, T(), T(), T(), T(), T(), T(), T(), (T)1, (T)-1, T(), T(), T(),
     T()},
    /*14*/
    {T(), T(), T(), T(), T(), T(), T(), T(), (T)1, (T)-1, T(), T(), T(), T(), T(), T(), T(), (T)1,
     (T)-1},
    /*15*/
    {T(), T(), T(), T(), T(), T(), (T)1, (T)-1, T(), T(), T(), T(), T(), T(), T(), (T)1, (T)-1, T(),
     T()},
    /*16*/
    {T(), T(), T(), T(), (T)-1, (T)-1, (T)1, (T)1, T(), T(), T(), T(), T(), (T)1, (T)1, (T)-1,
     (T)-1, T(), T()},
    /*17*/
    {T(), T(), T(), T(), (T)1, (T)-1, T(), T(), (T)-1, (T)-1, T(), T(), T(), (T)-1, (T)1, T(), T(),
     (T)1, (T)1},
    /*18*/
    {T(), T(), T(), T(), T(), T(), (T)-1, (T)1, (T)1, (T)-1, T(), T(), T(), T(), T(), (T)1, (T)-1,
     (T)-1, (T)1}};

template <typename T>
const T MRTD3Q19DescriptorBase<T>::invM[BaseDescriptor::q][BaseDescriptor::q] = {
    /*0*/ {
        (T)1 / (T)19, -(T)5 / (T)399, (T)1 / (T)21, T(), T(), T(), T(), T(), T(), T(), T(), T(),
        T(), T(), T(), T(), T(), T(), T()},

    /*1*/
    {(T)1 / (T)19, -(T)11 / (T)2394, -(T)1 / (T)63, -(T)1 / (T)10, (T)1 / (T)10, T(), T(), T(), T(),
     (T)1 / (T)18, -(T)1 / (T)18, T(), T(), T(), T(), T(), T(), T(), T()},

    /*2*/
    {(T)1 / (T)19, -(T)11 / (T)2394, -(T)1 / (T)63, T(), T(), -(T)1 / (T)10, (T)1 / (T)10, T(), T(),
     -(T)1 / (T)36, (T)1 / (T)36, (T)1 / (T)12, -(T)1 / (T)12, T(), T(), T(), T(), T(), T()},

    /*3*/
    {(T)1 / (T)19, -(T)11 / (T)2394, -(T)1 / (T)63, T(), T(), T(), T(), -(T)1 / (T)10, (T)1 / (T)10,
     -(T)1 / (T)36, (T)1 / (T)36, -(T)1 / (T)12, (T)1 / (T)12, T(), T(), T(), T(), T(), T()},

    /*4*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, -(T)1 / (T)10, -(T)1 / (T)40, -(T)1 / (T)10,
     -(T)1 / (T)40, T(), T(), (T)1 / (T)36, (T)1 / (T)72, (T)1 / (T)12, (T)1 / (T)24, (T)1 / (T)4,
     T(), T(), -(T)1 / (T)8, (T)1 / (T)8, T()},

    /*5*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, -(T)1 / (T)10, -(T)1 / (T)40, (T)1 / (T)10,
     (T)1 / (T)40, T(), T(), (T)1 / (T)36, (T)1 / (T)72, (T)1 / (T)12, (T)1 / (T)24, -(T)1 / (T)4,
     T(), T(), -(T)1 / (T)8, -(T)1 / (T)8, T()},

    /*6*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, -(T)1 / (T)10, -(T)1 / (T)40, T(), T(),
     -(T)1 / (T)10, -(T)1 / (T)40, (T)1 / (T)36, (T)1 / (T)72, -(T)1 / (T)12, -(T)1 / (T)24, T(),
     T(), (T)1 / (T)4, (T)1 / (T)8, T(), -(T)1 / (T)8},

    /*7*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, -(T)1 / (T)10, -(T)1 / (T)40, T(), T(),
     (T)1 / (T)10, (T)1 / (T)40, (T)1 / (T)36, (T)1 / (T)72, -(T)1 / (T)12, -(T)1 / (T)24, T(), T(),
     -(T)1 / (T)4, (T)1 / (T)8, T(), (T)1 / (T)8},

    /*8*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, T(), T(), -(T)1 / (T)10, -(T)1 / (T)40,
     -(T)1 / (T)10, -(T)1 / (T)40, -(T)1 / (T)18, -(T)1 / (T)36, T(), T(), T(), (T)1 / (T)4, T(),
     T(), -(T)1 / (T)8, (T)1 / (T)8},

    /*9*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, T(), T(), -(T)1 / (T)10, -(T)1 / (T)40,
     (T)1 / (T)10, (T)1 / (T)40, -(T)1 / (T)18, -(T)1 / (T)36, T(), T(), T(), -(T)1 / (T)4, T(),
     T(), -(T)1 / (T)8, -(T)1 / (T)8},

    /*10*/
    {(T)1 / (T)19, -(T)11 / (T)2394, -(T)1 / (T)63, (T)1 / (T)10, -(T)1 / (T)10, T(), T(), T(), T(),
     (T)1 / (T)18, -(T)1 / (T)18, T(), T(), T(), T(), T(), T(), T(), T()},

    /*11*/
    {(T)1 / (T)19, -(T)11 / (T)2394, -(T)1 / (T)63, T(), T(), (T)1 / (T)10, -(T)1 / (T)10, T(), T(),
     -(T)1 / (T)36, (T)1 / (T)36, (T)1 / (T)12, -(T)1 / (T)12, T(), T(), T(), T(), T(), T()},

    /*12*/
    {(T)1 / (T)19, -(T)11 / (T)2394, -(T)1 / (T)63, T(), T(), T(), T(), (T)1 / (T)10, -(T)1 / (T)10,
     -(T)1 / (T)36, (T)1 / (T)36, -(T)1 / (T)12, (T)1 / (T)12, T(), T(), T(), T(), T(), T()},

    /*13*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, (T)1 / (T)10, (T)1 / (T)40, (T)1 / (T)10,
     (T)1 / (T)40, T(), T(), (T)1 / (T)36, (T)1 / (T)72, (T)1 / (T)12, (T)1 / (T)24, (T)1 / (T)4,
     T(), T(), (T)1 / (T)8, -(T)1 / (T)8, T()},

    /*14*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, (T)1 / (T)10, (T)1 / (T)40, -(T)1 / (T)10,
     -(T)1 / (T)40, T(), T(), (T)1 / (T)36, (T)1 / (T)72, (T)1 / (T)12, (T)1 / (T)24, -(T)1 / (T)4,
     T(), T(), (T)1 / (T)8, (T)1 / (T)8, T()},

    /*15*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, (T)1 / (T)10, (T)1 / (T)40, T(), T(),
     (T)1 / (T)10, (T)1 / (T)40, (T)1 / (T)36, (T)1 / (T)72, -(T)1 / (T)12, -(T)1 / (T)24, T(), T(),
     (T)1 / (T)4, -(T)1 / (T)8, T(), (T)1 / (T)8},

    /*16*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, (T)1 / (T)10, (T)1 / (T)40, T(), T(),
     -(T)1 / (T)10, -(T)1 / (T)40, (T)1 / (T)36, (T)1 / (T)72, -(T)1 / (T)12, -(T)1 / (T)24, T(),
     T(), -(T)1 / (T)4, -(T)1 / (T)8, T(), -(T)1 / (T)8},

    /*17*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, T(), T(), (T)1 / (T)10, (T)1 / (T)40,
     (T)1 / (T)10, (T)1 / (T)40, -(T)1 / (T)18, -(T)1 / (T)36, T(), T(), T(), (T)1 / (T)4, T(), T(),
     (T)1 / (T)8, -(T)1 / (T)8},

    /*18*/
    {(T)1 / (T)19, (T)4 / (T)1197, (T)1 / (T)252, T(), T(), (T)1 / (T)10, (T)1 / (T)40,
     -(T)1 / (T)10, -(T)1 / (T)40, -(T)1 / (T)18, -(T)1 / (T)36, T(), T(), T(), -(T)1 / (T)4, T(),
     T(), (T)1 / (T)8, (T)1 / (T)8}};

template <typename T>
const T MRTD3Q19DescriptorBase<T>::S[BaseDescriptor::q] = {
    /*s0*/ (T)1,  // rho (conserved)
    /*s1*/ (T)1.19,
    /*s2*/ (T)1.4,
    /*s3*/ (T)1,    // rho*ux (conserved)
    /*s4*/ (T)1.2,  // <=> q index
    /*s5*/ (T)1,    // rho*uy (conserved)
    /*s6*/ (T)1.2,  // = s4 <=> q index
    /*s7*/ (T)1,    // rho*uz (conserved)
    /*s8*/ (T)1.2,  // = s4 <=> q index
    /*s9*/ T(),     // should be equal to s13, used to define nu
    /*s10*/ (T)1.4,
    /*s11*/ T(),  // = s9,
    /*s12*/ (T)1.4,
    /*s13*/ T(),  // should be equal to s9, used to define nu
    /*s14*/ T(),  // = s13,
    /*s15*/ T(),  // = s13,
    /*s16*/ (T)1.98,
    /*s17*/ (T)1.98,  // = s16,
    /*s18*/ (T)1.98   // = s16,
};

template <typename T>
const int MRTD3Q19DescriptorBase<T>::momentumIndexes[MRTD3Q19DescriptorBase<T>::jIndexes] = {
    3, 5, 7};

template <typename T>
const int MRTD3Q19DescriptorBase<T>::shearViscIndexes[MRTD3Q19DescriptorBase<T>::shearIndexes] = {
    9, 11, 13, 14, 15};

template <typename T>
const int MRTD3Q19DescriptorBase<T>::qViscIndexes[MRTD3Q19DescriptorBase<T>::qIndexes] = {4, 6, 8};

template <typename T>
const char MRTD3Q19Descriptor<T>::name[] = "MRTD3Q19";

// EXTERNAL FORCE MRT D3Q19 ////////////////////////////////////////////////////////////

template <typename T>
const char ForcedMRTD3Q19Descriptor<T>::name[] = "ForcedMRTD3Q19";

// External Force, density and momentum for ShanChen MRT 3D model
template <typename T>
const char ForcedMRTShanChenD3Q19Descriptor<T>::name[] = "ForcedShanChenMRTD3Q19";

// EXTERNAL VELOCITY MRT D3Q19 ////////////////////////////////////////////////////////////

template <typename T>
const char ExternalVelocityMRTD3Q19Descriptor<T>::name[] = "ExternalVelocityMRTD3Q19";

// HUI ABSORBING LAYER MRT D3Q19 ////////////////////////////////////////////////////////////

template <typename T>
const char AbsorbingWaveMRTD3Q19Descriptor<T>::name[] = "AbsorbingWaveMRTD3Q19";

// D3Q27 MRT
template <typename T>
const T MRTD3Q27DescriptorBase<T>::M[BaseDescriptor::q][BaseDescriptor::q] = {
    {(T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1,
     (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1},
    {T(),  -(T)1, T(), T(),  -(T)1, -(T)1, -(T)1, -(T)1, T(), T(),  -(T)1, -(T)1, -(T)1, -(T)1,
     (T)1, T(),   T(), (T)1, (T)1,  (T)1,  (T)1,  T(),   T(), (T)1, (T)1,  (T)1,  (T)1},
    {T(), T(),  -(T)1, T(),  -(T)1, (T)1, T(), T(),  -(T)1, -(T)1, -(T)1, -(T)1, (T)1, (T)1,
     T(), (T)1, T(),   (T)1, -(T)1, T(),  T(), (T)1, (T)1,  (T)1,  (T)1,  -(T)1, -(T)1},
    {T(), T(), T(),  -(T)1, T(), T(),  -(T)1, (T)1, -(T)1, (T)1, -(T)1, (T)1, -(T)1, (T)1,
     T(), T(), (T)1, T(),   T(), (T)1, -(T)1, (T)1, -(T)1, (T)1, -(T)1, (T)1, -(T)1},
    {-(T)2, -(T)1, -(T)1, -(T)1, T(), T(), T(), T(), T(), T(),  (T)1, (T)1, (T)1, (T)1,
     -(T)1, -(T)1, -(T)1, T(),   T(), T(), T(), T(), T(), (T)1, (T)1, (T)1, (T)1},
    {T(),  (T)2,  -(T)1, -(T)1, (T)1, (T)1, (T)1, (T)1,  -(T)2, -(T)2, T(), T(), T(), T(),
     (T)2, -(T)1, -(T)1, (T)1,  (T)1, (T)1, (T)1, -(T)2, -(T)2, T(),   T(), T(), T()},
    {T(), T(),  (T)1,  -(T)1, (T)1, (T)1,  -(T)1, -(T)1, T(), T(), T(), T(), T(), T(),
     T(), (T)1, -(T)1, (T)1,  (T)1, -(T)1, -(T)1, T(),   T(), T(), T(), T(), T()},
    {T(), T(), T(), T(),  (T)1,  -(T)1, T(), T(), T(), T(),  (T)1, (T)1,  -(T)1, -(T)1,
     T(), T(), T(), (T)1, -(T)1, T(),   T(), T(), T(), (T)1, (T)1, -(T)1, -(T)1},
    {T(), T(), T(), T(), T(), T(), T(), T(),  (T)1,  -(T)1, (T)1,  -(T)1, -(T)1, (T)1,
     T(), T(), T(), T(), T(), T(), T(), (T)1, -(T)1, (T)1,  -(T)1, -(T)1, (T)1},
    {T(), T(), T(), T(), T(), T(),  (T)1,  -(T)1, T(), T(),  (T)1,  -(T)1, (T)1, -(T)1,
     T(), T(), T(), T(), T(), (T)1, -(T)1, T(),   T(), (T)1, -(T)1, (T)1,  -(T)1},
    {T(),   (T)4, T(), T(),   (T)1,  (T)1,  (T)1,  (T)1, T(), T(),  -(T)2, -(T)2, -(T)2, -(T)2,
     -(T)4, T(),  T(), -(T)1, -(T)1, -(T)1, -(T)1, T(),  T(), (T)2, (T)2,  (T)2,  (T)2},
    {T(), T(),   (T)4, T(),   (T)1, -(T)1, T(), T(),   (T)1,  (T)1, -(T)2, -(T)2, (T)2, (T)2,
     T(), -(T)4, T(),  -(T)1, (T)1, T(),   T(), -(T)1, -(T)1, (T)2, (T)2,  -(T)2, -(T)2},
    {T(), T(), T(),   (T)4, T(), T(),   (T)1, -(T)1, (T)1, -(T)1, -(T)2, (T)2, -(T)2, (T)2,
     T(), T(), -(T)4, T(),  T(), -(T)1, (T)1, -(T)1, (T)1, (T)2,  -(T)2, (T)2, -(T)2},
    {T(),  -(T)4, T(), T(),   (T)2,  (T)2,  (T)2,  (T)2, T(), T(),  -(T)1, -(T)1, -(T)1, -(T)1,
     (T)4, T(),   T(), -(T)2, -(T)2, -(T)2, -(T)2, T(),  T(), (T)1, (T)1,  (T)1,  (T)1},
    {T(), T(),  -(T)4, T(),   (T)2, -(T)2, T(), T(),   (T)2,  (T)2, -(T)1, -(T)1, (T)1, (T)1,
     T(), (T)4, T(),   -(T)2, (T)2, T(),   T(), -(T)2, -(T)2, (T)1, (T)1,  -(T)1, -(T)1},
    {T(), T(), T(),  -(T)4, T(), T(),   (T)2, -(T)2, (T)2, -(T)2, -(T)1, (T)1, -(T)1, (T)1,
     T(), T(), (T)4, T(),   T(), -(T)2, (T)2, -(T)2, (T)2, (T)1,  -(T)1, (T)1, -(T)1},
    {(T)4, T(), T(), T(),   -(T)1, -(T)1, -(T)1, -(T)1, -(T)1, -(T)1, (T)1, (T)1, (T)1, (T)1,
     T(),  T(), T(), -(T)1, -(T)1, -(T)1, -(T)1, -(T)1, -(T)1, (T)1,  (T)1, (T)1, (T)1},
    {-(T)8, (T)4, (T)4, (T)4,  -(T)2, -(T)2, -(T)2, -(T)2, -(T)2, -(T)2, (T)1, (T)1, (T)1, (T)1,
     (T)4,  (T)4, (T)4, -(T)2, -(T)2, -(T)2, -(T)2, -(T)2, -(T)2, (T)1,  (T)1, (T)1, (T)1},
    {T(),   -(T)4, (T)2, (T)2, (T)1, (T)1, (T)1, (T)1,  -(T)2, -(T)2, T(), T(), T(), T(),
     -(T)4, (T)2,  (T)2, (T)1, (T)1, (T)1, (T)1, -(T)2, -(T)2, T(),   T(), T(), T()},
    {T(), T(),   -(T)2, (T)2, (T)1, (T)1,  -(T)1, -(T)1, T(), T(), T(), T(), T(), T(),
     T(), -(T)2, (T)2,  (T)1, (T)1, -(T)1, -(T)1, T(),   T(), T(), T(), T(), T()},
    {T(), T(), T(), T(),   -(T)2, (T)2, T(), T(), T(), T(),  (T)1, (T)1,  -(T)1, -(T)1,
     T(), T(), T(), -(T)2, (T)2,  T(),  T(), T(), T(), (T)1, (T)1, -(T)1, -(T)1},
    {T(), T(), T(), T(), T(), T(), T(), T(),   -(T)2, (T)2, (T)1,  -(T)1, -(T)1, (T)1,
     T(), T(), T(), T(), T(), T(), T(), -(T)2, (T)2,  (T)1, -(T)1, -(T)1, (T)1},
    {T(), T(), T(), T(), T(), T(),   -(T)2, (T)2, T(), T(),  (T)1,  -(T)1, (T)1, -(T)1,
     T(), T(), T(), T(), T(), -(T)2, (T)2,  T(),  T(), (T)1, -(T)1, (T)1,  -(T)1},
    {T(), T(), T(), T(),  -(T)1, -(T)1, (T)1,  (T)1, T(), T(), T(), T(), T(), T(),
     T(), T(), T(), (T)1, (T)1,  -(T)1, -(T)1, T(),  T(), T(), T(), T(), T()},
    {T(), T(), T(), T(),   (T)1, -(T)1, T(), T(),  -(T)1, -(T)1, T(), T(), T(), T(),
     T(), T(), T(), -(T)1, (T)1, T(),   T(), (T)1, (T)1,  T(),   T(), T(), T()},
    {T(), T(), T(), T(), T(), T(),  -(T)1, (T)1,  (T)1, -(T)1, T(), T(), T(), T(),
     T(), T(), T(), T(), T(), (T)1, -(T)1, -(T)1, (T)1, T(),   T(), T(), T()},
    {T(), T(), T(), T(), T(), T(), T(), T(), T(), T(),  -(T)1, (T)1,  (T)1, -(T)1,
     T(), T(), T(), T(), T(), T(), T(), T(), T(), (T)1, -(T)1, -(T)1, (T)1}};

template <typename T>
const T MRTD3Q27DescriptorBase<T>::invM[BaseDescriptor::q][BaseDescriptor::q] = {
    {(T)1 / (T)27, T(), T(), T(), -(T)1 / (T)9, T(), T(), T(),         T(),
     T(),          T(), T(), T(), T(),          T(), T(), (T)1 / (T)9, -(T)1 / (T)27,
     T(),          T(), T(), T(), T(),          T(), T(), T(),         T()},
    {(T)1 / (T)27,  -(T)1 / (T)18, T(), T(), -(T)1 / (T)18, (T)1 / (T)18, T(), T(), T(),
     T(),           (T)1 / (T)18,  T(), T(), -(T)1 / (T)18, T(),          T(), T(), (T)1 / (T)54,
     -(T)1 / (T)18, T(),           T(), T(), T(),           T(),          T(), T(), T()},
    {(T)1 / (T)27,  T(),          -(T)1 / (T)18, T(),          -(T)1 / (T)18,
     -(T)1 / (T)36, (T)1 / (T)12, T(),           T(),          T(),
     T(),           (T)1 / (T)18, T(),           T(),          -(T)1 / (T)18,
     T(),           T(),          (T)1 / (T)54,  (T)1 / (T)36, -(T)1 / (T)12,
     T(),           T(),          T(),           T(),          T(),
     T(),           T()},
    {(T)1 / (T)27,  T(),          T(), -(T)1 / (T)18, -(T)1 / (T)18, -(T)1 / (T)36,
     -(T)1 / (T)12, T(),          T(), T(),           T(),           T(),
     (T)1 / (T)18,  T(),          T(), -(T)1 / (T)18, T(),           (T)1 / (T)54,
     (T)1 / (T)36,  (T)1 / (T)12, T(), T(),           T(),           T(),
     T(),           T(),          T()},
    {(T)1 / (T)27, -(T)1 / (T)18, -(T)1 / (T)18,
     T(),          T(),           (T)1 / (T)36,
     (T)1 / (T)12, (T)1 / (T)12,  T(),
     T(),          (T)1 / (T)72,  (T)1 / (T)72,
     T(),          (T)1 / (T)36,  (T)1 / (T)36,
     T(),          -(T)1 / (T)36, -(T)1 / (T)108,
     (T)1 / (T)72, (T)1 / (T)24,  -(T)1 / (T)12,
     T(),          T(),           -(T)1 / (T)8,
     (T)1 / (T)8,  T(),           T()},
    {(T)1 / (T)27, -(T)1 / (T)18, (T)1 / (T)18,
     T(),          T(),           (T)1 / (T)36,
     (T)1 / (T)12, -(T)1 / (T)12, T(),
     T(),          (T)1 / (T)72,  -(T)1 / (T)72,
     T(),          (T)1 / (T)36,  -(T)1 / (T)36,
     T(),          -(T)1 / (T)36, -(T)1 / (T)108,
     (T)1 / (T)72, (T)1 / (T)24,  (T)1 / (T)12,
     T(),          T(),           -(T)1 / (T)8,
     -(T)1 / (T)8, T(),           T()},
    {(T)1 / (T)27,  -(T)1 / (T)18, T(), -(T)1 / (T)18, T(),           (T)1 / (T)36,
     -(T)1 / (T)12, T(),           T(), (T)1 / (T)12,  (T)1 / (T)72,  T(),
     (T)1 / (T)72,  (T)1 / (T)36,  T(), (T)1 / (T)36,  -(T)1 / (T)36, -(T)1 / (T)108,
     (T)1 / (T)72,  -(T)1 / (T)24, T(), T(),           -(T)1 / (T)12, (T)1 / (T)8,
     T(),           -(T)1 / (T)8,  T()},
    {(T)1 / (T)27,  -(T)1 / (T)18, T(), (T)1 / (T)18,  T(),           (T)1 / (T)36,
     -(T)1 / (T)12, T(),           T(), -(T)1 / (T)12, (T)1 / (T)72,  T(),
     -(T)1 / (T)72, (T)1 / (T)36,  T(), -(T)1 / (T)36, -(T)1 / (T)36, -(T)1 / (T)108,
     (T)1 / (T)72,  -(T)1 / (T)24, T(), T(),           (T)1 / (T)12,  (T)1 / (T)8,
     T(),           (T)1 / (T)8,   T()},
    {(T)1 / (T)27,  T(),           -(T)1 / (T)18,
     -(T)1 / (T)18, T(),           -(T)1 / (T)18,
     T(),           T(),           (T)1 / (T)12,
     T(),           T(),           (T)1 / (T)72,
     (T)1 / (T)72,  T(),           (T)1 / (T)36,
     (T)1 / (T)36,  -(T)1 / (T)36, -(T)1 / (T)108,
     -(T)1 / (T)36, T(),           T(),
     -(T)1 / (T)12, T(),           T(),
     -(T)1 / (T)8,  (T)1 / (T)8,   T()},
    {(T)1 / (T)27,  T(),           -(T)1 / (T)18,
     (T)1 / (T)18,  T(),           -(T)1 / (T)18,
     T(),           T(),           -(T)1 / (T)12,
     T(),           T(),           (T)1 / (T)72,
     -(T)1 / (T)72, T(),           (T)1 / (T)36,
     -(T)1 / (T)36, -(T)1 / (T)36, -(T)1 / (T)108,
     -(T)1 / (T)36, T(),           T(),
     (T)1 / (T)12,  T(),           T(),
     -(T)1 / (T)8,  -(T)1 / (T)8,  T()},
    {(T)1 / (T)27,  -(T)1 / (T)18, -(T)1 / (T)18, -(T)1 / (T)18, (T)1 / (T)18,  T(),
     T(),           (T)1 / (T)12,  (T)1 / (T)12,  (T)1 / (T)12,  -(T)1 / (T)36, -(T)1 / (T)36,
     -(T)1 / (T)36, -(T)1 / (T)72, -(T)1 / (T)72, -(T)1 / (T)72, (T)1 / (T)36,  (T)1 / (T)216,
     T(),           T(),           (T)1 / (T)24,  (T)1 / (T)24,  (T)1 / (T)24,  T(),
     T(),           T(),           -(T)1 / (T)8},
    {(T)1 / (T)27, -(T)1 / (T)18, -(T)1 / (T)18, (T)1 / (T)18,  (T)1 / (T)18,  T(),
     T(),          (T)1 / (T)12,  -(T)1 / (T)12, -(T)1 / (T)12, -(T)1 / (T)36, -(T)1 / (T)36,
     (T)1 / (T)36, -(T)1 / (T)72, -(T)1 / (T)72, (T)1 / (T)72,  (T)1 / (T)36,  (T)1 / (T)216,
     T(),          T(),           (T)1 / (T)24,  -(T)1 / (T)24, -(T)1 / (T)24, T(),
     T(),          T(),           (T)1 / (T)8},
    {(T)1 / (T)27,  -(T)1 / (T)18, (T)1 / (T)18,  -(T)1 / (T)18, (T)1 / (T)18,  T(),
     T(),           -(T)1 / (T)12, -(T)1 / (T)12, (T)1 / (T)12,  -(T)1 / (T)36, (T)1 / (T)36,
     -(T)1 / (T)36, -(T)1 / (T)72, (T)1 / (T)72,  -(T)1 / (T)72, (T)1 / (T)36,  (T)1 / (T)216,
     T(),           T(),           -(T)1 / (T)24, -(T)1 / (T)24, (T)1 / (T)24,  T(),
     T(),           T(),           (T)1 / (T)8},
    {(T)1 / (T)27, -(T)1 / (T)18, (T)1 / (T)18,  (T)1 / (T)18,  (T)1 / (T)18,  T(),
     T(),          -(T)1 / (T)12, (T)1 / (T)12,  -(T)1 / (T)12, -(T)1 / (T)36, (T)1 / (T)36,
     (T)1 / (T)36, -(T)1 / (T)72, (T)1 / (T)72,  (T)1 / (T)72,  (T)1 / (T)36,  (T)1 / (T)216,
     T(),          T(),           -(T)1 / (T)24, (T)1 / (T)24,  -(T)1 / (T)24, T(),
     T(),          T(),           -(T)1 / (T)8},
    {(T)1 / (T)27,  (T)1 / (T)18,  T(), T(), -(T)1 / (T)18, (T)1 / (T)18, T(), T(), T(),
     T(),           -(T)1 / (T)18, T(), T(), (T)1 / (T)18,  T(),          T(), T(), (T)1 / (T)54,
     -(T)1 / (T)18, T(),           T(), T(), T(),           T(),          T(), T(), T()},
    {(T)1 / (T)27, T(),           (T)1 / (T)18, T(), -(T)1 / (T)18, -(T)1 / (T)36,
     (T)1 / (T)12, T(),           T(),          T(), T(),           -(T)1 / (T)18,
     T(),          T(),           (T)1 / (T)18, T(), T(),           (T)1 / (T)54,
     (T)1 / (T)36, -(T)1 / (T)12, T(),          T(), T(),           T(),
     T(),          T(),           T()},
    {(T)1 / (T)27,  T(),          T(), (T)1 / (T)18, -(T)1 / (T)18, -(T)1 / (T)36,
     -(T)1 / (T)12, T(),          T(), T(),          T(),           T(),
     -(T)1 / (T)18, T(),          T(), (T)1 / (T)18, T(),           (T)1 / (T)54,
     (T)1 / (T)36,  (T)1 / (T)12, T(), T(),          T(),           T(),
     T(),           T(),          T()},
    {(T)1 / (T)27, (T)1 / (T)18,  (T)1 / (T)18,
     T(),          T(),           (T)1 / (T)36,
     (T)1 / (T)12, (T)1 / (T)12,  T(),
     T(),          -(T)1 / (T)72, -(T)1 / (T)72,
     T(),          -(T)1 / (T)36, -(T)1 / (T)36,
     T(),          -(T)1 / (T)36, -(T)1 / (T)108,
     (T)1 / (T)72, (T)1 / (T)24,  -(T)1 / (T)12,
     T(),          T(),           (T)1 / (T)8,
     -(T)1 / (T)8, T(),           T()},
    {(T)1 / (T)27, (T)1 / (T)18,  -(T)1 / (T)18,
     T(),          T(),           (T)1 / (T)36,
     (T)1 / (T)12, -(T)1 / (T)12, T(),
     T(),          -(T)1 / (T)72, (T)1 / (T)72,
     T(),          -(T)1 / (T)36, (T)1 / (T)36,
     T(),          -(T)1 / (T)36, -(T)1 / (T)108,
     (T)1 / (T)72, (T)1 / (T)24,  (T)1 / (T)12,
     T(),          T(),           (T)1 / (T)8,
     (T)1 / (T)8,  T(),           T()},
    {(T)1 / (T)27,  (T)1 / (T)18,  T(), (T)1 / (T)18,  T(),           (T)1 / (T)36,
     -(T)1 / (T)12, T(),           T(), (T)1 / (T)12,  -(T)1 / (T)72, T(),
     -(T)1 / (T)72, -(T)1 / (T)36, T(), -(T)1 / (T)36, -(T)1 / (T)36, -(T)1 / (T)108,
     (T)1 / (T)72,  -(T)1 / (T)24, T(), T(),           -(T)1 / (T)12, -(T)1 / (T)8,
     T(),           (T)1 / (T)8,   T()},
    {(T)1 / (T)27,  (T)1 / (T)18,  T(), -(T)1 / (T)18, T(),           (T)1 / (T)36,
     -(T)1 / (T)12, T(),           T(), -(T)1 / (T)12, -(T)1 / (T)72, T(),
     (T)1 / (T)72,  -(T)1 / (T)36, T(), (T)1 / (T)36,  -(T)1 / (T)36, -(T)1 / (T)108,
     (T)1 / (T)72,  -(T)1 / (T)24, T(), T(),           (T)1 / (T)12,  -(T)1 / (T)8,
     T(),           -(T)1 / (T)8,  T()},
    {(T)1 / (T)27,  T(),           (T)1 / (T)18,
     (T)1 / (T)18,  T(),           -(T)1 / (T)18,
     T(),           T(),           (T)1 / (T)12,
     T(),           T(),           -(T)1 / (T)72,
     -(T)1 / (T)72, T(),           -(T)1 / (T)36,
     -(T)1 / (T)36, -(T)1 / (T)36, -(T)1 / (T)108,
     -(T)1 / (T)36, T(),           T(),
     -(T)1 / (T)12, T(),           T(),
     (T)1 / (T)8,   -(T)1 / (T)8,  T()},
    {(T)1 / (T)27,  T(),           (T)1 / (T)18,
     -(T)1 / (T)18, T(),           -(T)1 / (T)18,
     T(),           T(),           -(T)1 / (T)12,
     T(),           T(),           -(T)1 / (T)72,
     (T)1 / (T)72,  T(),           -(T)1 / (T)36,
     (T)1 / (T)36,  -(T)1 / (T)36, -(T)1 / (T)108,
     -(T)1 / (T)36, T(),           T(),
     (T)1 / (T)12,  T(),           T(),
     (T)1 / (T)8,   (T)1 / (T)8,   T()},
    {(T)1 / (T)27, (T)1 / (T)18, (T)1 / (T)18, (T)1 / (T)18, (T)1 / (T)18, T(),
     T(),          (T)1 / (T)12, (T)1 / (T)12, (T)1 / (T)12, (T)1 / (T)36, (T)1 / (T)36,
     (T)1 / (T)36, (T)1 / (T)72, (T)1 / (T)72, (T)1 / (T)72, (T)1 / (T)36, (T)1 / (T)216,
     T(),          T(),          (T)1 / (T)24, (T)1 / (T)24, (T)1 / (T)24, T(),
     T(),          T(),          (T)1 / (T)8},
    {(T)1 / (T)27,  (T)1 / (T)18, (T)1 / (T)18,  -(T)1 / (T)18, (T)1 / (T)18,  T(),
     T(),           (T)1 / (T)12, -(T)1 / (T)12, -(T)1 / (T)12, (T)1 / (T)36,  (T)1 / (T)36,
     -(T)1 / (T)36, (T)1 / (T)72, (T)1 / (T)72,  -(T)1 / (T)72, (T)1 / (T)36,  (T)1 / (T)216,
     T(),           T(),          (T)1 / (T)24,  -(T)1 / (T)24, -(T)1 / (T)24, T(),
     T(),           T(),          -(T)1 / (T)8},
    {(T)1 / (T)27, (T)1 / (T)18,  -(T)1 / (T)18, (T)1 / (T)18,  (T)1 / (T)18, T(),
     T(),          -(T)1 / (T)12, -(T)1 / (T)12, (T)1 / (T)12,  (T)1 / (T)36, -(T)1 / (T)36,
     (T)1 / (T)36, (T)1 / (T)72,  -(T)1 / (T)72, (T)1 / (T)72,  (T)1 / (T)36, (T)1 / (T)216,
     T(),          T(),           -(T)1 / (T)24, -(T)1 / (T)24, (T)1 / (T)24, T(),
     T(),          T(),           -(T)1 / (T)8},
    {(T)1 / (T)27,  (T)1 / (T)18,  -(T)1 / (T)18, -(T)1 / (T)18, (T)1 / (T)18,  T(),
     T(),           -(T)1 / (T)12, (T)1 / (T)12,  -(T)1 / (T)12, (T)1 / (T)36,  -(T)1 / (T)36,
     -(T)1 / (T)36, (T)1 / (T)72,  -(T)1 / (T)72, -(T)1 / (T)72, (T)1 / (T)36,  (T)1 / (T)216,
     T(),           T(),           -(T)1 / (T)24, (T)1 / (T)24,  -(T)1 / (T)24, T(),
     T(),           T(),           (T)1 / (T)8}};

template <typename T>
const int MRTD3Q27DescriptorBase<T>::momentumIndexes[MRTD3Q27DescriptorBase<T>::jIndexes] = {
    1, 2, 3};

template <typename T>
const int MRTD3Q27DescriptorBase<T>::shearViscIndexes[MRTD3Q27DescriptorBase<T>::shearIndexes] = {
    5, 6, 7, 8, 9};

template <typename T>
const T MRTD3Q27DescriptorBase<T>::S[BaseDescriptor::q] = {
    /*s0*/ (T)1,     // rho (conserved)
    /*s1*/ (T)1.1,   // unknown value. set to 1.1 arbitrarily
    /*s2*/ (T)1,     // j_x
    /*s3*/ (T)1,     // j_y
    /*s4*/ (T)1,     // j_z
    /*s5*/ T(),      // <=> nu index
    /*s6*/ T(),      // = s5
    /*s7*/ T(),      // = s5
    /*s8*/ T(),      // = s5
    /*s9*/ T(),      // = s5
    /*s10*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s11*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s12*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s13*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s14*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s15*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s16*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s17*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s18*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s19*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s20*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s21*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s22*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s23*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s24*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s25*/ (T)1.1,  // unknown value. set to 1.1 arbitrarily
    /*s26*/ (T)1.1   // unknown value. set to 1.1 arbitrarily
};

template <typename T>
const char MRTD3Q27Descriptor<T>::name[] = "MRTD3Q27";

}  // namespace descriptors

}  // namespace plb

#endif
