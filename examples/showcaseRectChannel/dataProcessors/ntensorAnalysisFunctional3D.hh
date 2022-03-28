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
 * Helper functions for domain initialization -- header file.
 */
#ifndef NTENSOR_ANALYSIS_FUNCTIONAL_3D_HH
#define NTENSOR_ANALYSIS_FUNCTIONAL_3D_HH

#include <cmath>
#include <limits>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "dataProcessors/ntensorAnalysisFunctional3D.h"
#include "finiteDifference/fdStencils1D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

template <typename T1, typename T2>
void CopyConvertNTensorFunctional3D<T1, T2>::process(
    Box3D domain, NTensorField3D<T1> &field1, NTensorField3D<T2> &field2)
{
    PLB_PRECONDITION(field1.getNdim() == field2.getNdim());
    plint ndim = field1.getNdim();
    Dot3D offset = computeRelativeDisplacement(field1, field2);

    // Improve computational speed if the field is scalar-valued.
    if (ndim == 1) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    for (int iDim = 0; iDim < ndim; ++iDim) {
                        *field2.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                            (T2)*field1.get(iX, iY, iZ);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else
    {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    for (int iDim = 0; iDim < ndim; ++iDim) {
                        field2.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iDim] =
                            (T2)field1.get(iX, iY, iZ)[iDim];
                    }
                }
            }
        }
    }
}

template <typename T1, typename T2>
CopyConvertNTensorFunctional3D<T1, T2> *CopyConvertNTensorFunctional3D<T1, T2>::clone() const
{
    return new CopyConvertNTensorFunctional3D<T1, T2>(*this);
}

template <typename T1, typename T2>
void CopyConvertNTensorFunctional3D<T1, T2>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T1, typename T2>
BlockDomain::DomainT CopyConvertNTensorFunctional3D<T1, T2>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ConvertNTensorToScalarFunctional3D ************************************ */

template <typename T>
void ConvertNTensorToScalarFunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 2);
    NTensorField3D<T> &nTensor = *dynamic_cast<NTensorField3D<T> *>(fields[0]);
    ScalarField3D<T> &scalar = *dynamic_cast<ScalarField3D<T> *>(fields[1]);
    PLB_ASSERT(nTensor.getNdim() == 1);

    Dot3D offset = computeRelativeDisplacement(nTensor, scalar);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                scalar.get(iX + offset.x, iY + offset.y, iZ + offset.z) = *nTensor.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T>
ConvertNTensorToScalarFunctional3D<T> *ConvertNTensorToScalarFunctional3D<T>::clone() const
{
    return new ConvertNTensorToScalarFunctional3D<T>(*this);
}

template <typename T>
void ConvertNTensorToScalarFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // nTensor
    modified[1] = modif::staticVariables;  // scalar
}

template <typename T>
BlockDomain::DomainT ConvertNTensorToScalarFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ConvertScalarToNTensorFunctional3D ************************************ */

template <typename T>
void ConvertScalarToNTensorFunctional3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> fields)
{
    PLB_PRECONDITION(fields.size() == 2);
    ScalarField3D<T> &scalar = *dynamic_cast<ScalarField3D<T> *>(fields[0]);
    NTensorField3D<T> &nTensor = *dynamic_cast<NTensorField3D<T> *>(fields[1]);
    PLB_ASSERT(nTensor.getNdim() == 1);

    Dot3D offset = computeRelativeDisplacement(scalar, nTensor);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                *nTensor.get(iX + offset.x, iY + offset.y, iZ + offset.z) = scalar.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T>
ConvertScalarToNTensorFunctional3D<T> *ConvertScalarToNTensorFunctional3D<T>::clone() const
{
    return new ConvertScalarToNTensorFunctional3D<T>(*this);
}

template <typename T>
void ConvertScalarToNTensorFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // scalar
    modified[1] = modif::staticVariables;  // nTensor
}

template <typename T>
BlockDomain::DomainT ConvertScalarToNTensorFunctional3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // NTENSOR_ANALYSIS_FUNCTIONAL_3D_HH
