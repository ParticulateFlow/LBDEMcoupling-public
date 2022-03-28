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
#ifndef NTENSOR_ANALYSIS_FUNCTIONAL_2D_HH
#define NTENSOR_ANALYSIS_FUNCTIONAL_2D_HH

#include <cmath>
#include <limits>

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "dataProcessors/ntensorAnalysisFunctional2D.h"
#include "finiteDifference/fdStencils1D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

template <typename T1, typename T2>
void CopyConvertNTensorFunctional2D<T1, T2>::process(
    Box2D domain, NTensorField2D<T1> &field1, NTensorField2D<T2> &field2)
{
    PLB_PRECONDITION(field1.getNdim() == field2.getNdim());
    plint ndim = field1.getNdim();
    Dot2D offset = computeRelativeDisplacement(field1, field2);

    // Improve computational speed if the field is scalar-valued.
    if (ndim == 1) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (int iDim = 0; iDim < ndim; ++iDim) {
                    *field2.get(iX + offset.x, iY + offset.y) = (T2)*field1.get(iX, iY);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else
    {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (int iDim = 0; iDim < ndim; ++iDim) {
                    field2.get(iX + offset.x, iY + offset.y)[iDim] = (T2)field1.get(iX, iY)[iDim];
                }
            }
        }
    }
}

template <typename T1, typename T2>
CopyConvertNTensorFunctional2D<T1, T2> *CopyConvertNTensorFunctional2D<T1, T2>::clone() const
{
    return new CopyConvertNTensorFunctional2D<T1, T2>(*this);
}

template <typename T1, typename T2>
void CopyConvertNTensorFunctional2D<T1, T2>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T1, typename T2>
BlockDomain::DomainT CopyConvertNTensorFunctional2D<T1, T2>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // NTENSOR_ANALYSIS_FUNCTIONAL_2D_HH
