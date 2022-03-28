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
 *  -- header file
 */
#ifndef MRT_LATTICES_H
#define MRT_LATTICES_H

#include <vector>

#include "core/globalDefs.h"
#include "latticeBoltzmann/externalFields.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "multiPhysics/shanChenLattices2D.h"
#include "multiPhysics/shanChenLattices3D.h"

namespace plb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
 */

/* Sébastien Leclaire contributed for the D3Q15 and D2Q27 lattices definitions */

namespace descriptors {

/// MRT D2Q9 lattice. The numbering follows the one in "Viscous flow computations
/// with the method of lattice Boltzmann equation", D. Yu, L.-S. Luo, W. Shi,
/// Progress in Aerospace Sciences 39, (2003), p. 329-367
template <typename T>
struct MRTD2Q9DescriptorBase : public D2Q9DescriptorBase<T> {
    typedef D2Q9DescriptorBase<T> BaseDescriptor;
    typedef MRTD2Q9DescriptorBase<T> SecondBaseDescriptor;
    enum { numPop = BaseDescriptor::q };

    static const T M[BaseDescriptor::q][BaseDescriptor::q];  // Matrix of base change between f and
                                                             // moments : moments=M.f
    static const T invM[BaseDescriptor::q]
                       [BaseDescriptor::q];  // inverse of base change matrix : f=invM.moments
    static const T S[BaseDescriptor::q];     // relaxation times
    enum { jIndexes = 2 };
    static const int momentumIndexes[jIndexes];  // relevant indexes of r. t. for shear viscosity
    enum { shearIndexes = 2 };
    static const int
        shearViscIndexes[shearIndexes];  // relevant indexes of r. t. for shear viscosity
    enum { qIndexes = 2 };
    static const int qViscIndexes[qIndexes];  // relevant indexes of r. t. for q
    static const int bulkViscIndex = 1;       // relevant index of r. t. for bulk viscosity
    static const int epsilonIndex = 2;        // relevant index of r. t. for epsilon
};

// Humières, D. (2002). Multiple–relaxation–time lattice Boltzmann models in three dimensions.
// Philosophical Transactions of the Royal Society of London A: Mathematical, Physical and
// Engineering Sciences, 360(1792), 437–451. Retrieved from
// http://rsta.royalsocietypublishing.org/content/360/1792/437.abstract
template <typename T>
struct MRTD3Q15DescriptorBase : public D3Q15DescriptorBase<T> {
    typedef D3Q15DescriptorBase<T> BaseDescriptor;
    typedef MRTD3Q15DescriptorBase<T> SecondBaseDescriptor;
    enum { numPop = BaseDescriptor::q };

    static const T M[BaseDescriptor::q][BaseDescriptor::q];  // Matrix of base change between f and
                                                             // moments : moments=M.f
    static const T invM[BaseDescriptor::q]
                       [BaseDescriptor::q];  // inverse of base change matrix : f=invM.moments
    static const T S[BaseDescriptor::q];     // relaxation times
    enum { jIndexes = 3 };
    static const int momentumIndexes[jIndexes];  // relevant indexes of r. t. for shear viscosity
    enum { shearIndexes = 5 };
    static const int
        shearViscIndexes[shearIndexes];  // relevant indexes of r. t. for shear viscosity
    static const int bulkViscIndex = 1;  // relevant index of r. t. for bulk viscosity
    enum { qIndexes = 3 };
    static const int
        qViscIndexes[qIndexes];  // relevant indexes of r. t. for q indices (s4,s6,s8 of the paper)
    static const int epsilonIndex =
        2;                         // relevant index of r. t. for epsilon (s2 in the original paper)
    static const int mIndex = 14;  // relevant index of r. t. for m (s14 in the original paper)
};

/// MRT D3Q19 lattice. The numbering follows the one in "Multiple-relaxation-
/// time lattice Boltzmann models in three dimensions", D. D'Humières,
/// I. Ginzburg, M. Krafzcyk, P. Lallemand, L.-S. Luo,
/// Phil. Trans. R. Soc. Lond. A (2002) 660, p. 437-451
template <typename T>
struct MRTD3Q19DescriptorBase : public D3Q19DescriptorBase<T> {
    typedef D3Q19DescriptorBase<T> BaseDescriptor;
    typedef MRTD3Q19DescriptorBase<T> SecondBaseDescriptor;
    enum { numPop = BaseDescriptor::q };

    static const T M[BaseDescriptor::q][BaseDescriptor::q];  // Matrix of base change between f and
                                                             // moments : moments=M.f
    static const T invM[BaseDescriptor::q]
                       [BaseDescriptor::q];  // inverse of base change matrix : f=invM.moments
    static const T S[BaseDescriptor::q];     // relaxation times
    enum { jIndexes = 3 };
    static const int momentumIndexes[jIndexes];  // relevant indexes of r. t. for shear viscosity
    enum { shearIndexes = 5 };
    static const int
        shearViscIndexes[shearIndexes];  // relevant indexes of r. t. for shear viscosity
    static const int bulkViscIndex = 1;  // relevant index of r. t. for bulk viscosity
    enum { qIndexes = 3 };
    static const int
        qViscIndexes[qIndexes];  // relevant indexes of r. t. for q indices (s4,s6,s8 of the paper)
    static const int epsilonIndex =
        2;  // relevant index of r. t. for epsilon (s2 in the original paper)
};

/// MRT D3Q27 lattice.
/// Dubois, F., & Lallemand, P. (2011). Quartic parameters for acoustic applications of lattice
/// Boltzmann scheme. Computers & Mathematics with Applications, 61(12), 3404–3416.
/// doi:10.1016/j.camwa.2011.01.011

/// Suga, K., Kuwata, Y., Takashima, K., & Chikasue, R. (2015). A D3Q27 multiple-relaxation-time
/// lattice Boltzmann method for turbulent flows. Computers & Mathematics with Applications, 69(6),
/// 518–529. doi:10.1016/j.camwa.2015.01.010
template <typename T>
struct MRTD3Q27DescriptorBase : public D3Q27DescriptorBase<T> {
    typedef D3Q27DescriptorBase<T> BaseDescriptor;
    typedef MRTD3Q27DescriptorBase<T> SecondBaseDescriptor;
    enum { numPop = BaseDescriptor::q };

    static const T M[BaseDescriptor::q][BaseDescriptor::q];  // Matrix of base change between f and
                                                             // moments : moments=M.f
    static const T invM[BaseDescriptor::q]
                       [BaseDescriptor::q];  // inverse of base change matrix : f=invM.moments
    static const T S[BaseDescriptor::q];     // relaxation times
    enum { jIndexes = 3 };
    static const int momentumIndexes[jIndexes];  // relevant indexes of r. t. for shear viscosity
    enum { shearIndexes = 5 };
    static const int
        shearViscIndexes[shearIndexes];  // relevant indexes of r. t. for shear viscosity
    static const int bulkViscIndex = 1;  // relevant index of r. t. for bulk viscosity
    enum { qIndexes = 3 };
    static const int
        qViscIndexes[qIndexes];  // relevant indexes of r. t. for q indices (s4,s6,s8 of the paper)
    static const int epsilonIndex =
        2;  // relevant index of r. t. for epsilon (s2 in the original paper)
};

template <typename T>
struct MRTD2Q9Descriptor : public MRTD2Q9DescriptorBase<T>, public NoExternalFieldBase {
    static const char name[];
};

template <typename T>
struct ForcedMRTD2Q9Descriptor : public MRTD2Q9DescriptorBase<T>, public Force2dDescriptorBase {
    static const char name[];
};

template <typename T>
struct ExternalVelocityMRTD2Q9Descriptor :
    public MRTD2Q9DescriptorBase<T>,
    public Velocity2dDescriptorBase {
    static const char name[];
};

template <typename T>
struct ForcedMRTShanChenD2Q9Descriptor :
    public MRTD2Q9DescriptorBase<T>,
    public ForcedShanChenExternalBase2D {
    static const char name[];
};

template <typename T>
struct MRTD3Q15Descriptor : public MRTD3Q15DescriptorBase<T>, public NoExternalFieldBase {
    static const char name[];
};

template <typename T>
struct MRTD3Q19Descriptor : public MRTD3Q19DescriptorBase<T>, public NoExternalFieldBase {
    static const char name[];
};

template <typename T>
struct ForcedMRTD3Q19Descriptor : public MRTD3Q19DescriptorBase<T>, public Force3dDescriptorBase {
    static const char name[];
};

template <typename T>
struct ForcedMRTShanChenD3Q19Descriptor :
    public MRTD3Q19DescriptorBase<T>,
    public ForcedShanChenExternalBase3D {
    static const char name[];
};

template <typename T>
struct ExternalVelocityMRTD3Q19Descriptor :
    public MRTD3Q19DescriptorBase<T>,
    public Velocity3dBase {
    static const char name[];
};

template <typename T>
struct AbsorbingWaveMRTD3Q19Descriptor :
    public MRTD3Q19DescriptorBase<T>,
    public AbsorbingWaveExternalField3dBase {
    static const char name[];
};

template <typename T>
struct MRTD3Q27Descriptor : public MRTD3Q27DescriptorBase<T>, public NoExternalFieldBase {
    static const char name[];
};

}  // namespace descriptors

}  // namespace plb

#endif
