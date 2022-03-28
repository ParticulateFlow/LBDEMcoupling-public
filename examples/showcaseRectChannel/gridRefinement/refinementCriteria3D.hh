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

#ifndef REFINEMENT_CRITERIA_3D_HH
#define REFINEMENT_CRITERIA_3D_HH

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/util.h"
#include "gridRefinement/refinementCriteria3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
ComputeRefinementRvalueFunctional3D<T, Descriptor>::ComputeRefinementRvalueFunctional3D(
    T knudsen_) :
    knudsen(knudsen_)
{ }

template <typename T, template <typename U> class Descriptor>
void ComputeRefinementRvalueFunctional3D<T, Descriptor>::process(
    Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<T> &scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> const &cell = lattice.get(iX, iY, iZ);
                Array<T, 3> j;
                T rhoBar;
                cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);
                T jSqr = normSqr(j);
                T C = (T)0;
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    T f = cell[iPop];
                    T feq = cell.getDynamics().computeEquilibrium(iPop, rhoBar, j, jSqr);
                    T fneq = f - feq;
                    C += std::fabs(
                        fneq
                        / fullF<T, Descriptor>(
                            feq, iPop));  // remember : f may be equal to f - t[i]
                }
                C /= ((T)Descriptor<T>::q * knudsen);
                if (C <= (T)0) {  // if C <= 0 then there is no need for refinment here
                    scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) = (T)0;
                } else {
                    T val = (T)std::log(C) / (T)std::log((T)2);
                    if (val < (T)0) {
                        val = (T)0;
                    }
                    scalarField.get(iX + offset.x, iY + offset.y, iZ + offset.z) = val;
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
ComputeRefinementRvalueFunctional3D<T, Descriptor>
    *ComputeRefinementRvalueFunctional3D<T, Descriptor>::clone() const
{
    return new ComputeRefinementRvalueFunctional3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void ComputeRefinementRvalueFunctional3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // Lattice.
    modified[1] = modif::staticVariables;  // Scalar Field.
}

template <typename T, template <typename U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > computeRvalues(
    MultiBlockLattice3D<T, Descriptor> &lattice, T knudsen, Box3D const &domain)
{
    std::unique_ptr<MultiScalarField3D<T> > rValues = generateMultiScalarField<T>(lattice, domain);
    applyProcessingFunctional(
        new ComputeRefinementRvalueFunctional3D<T, Descriptor>(knudsen), domain, lattice, *rValues);

    return rValues;
}

}  // namespace plb

#endif  // REFINEMENT_CRITERIA_3D_HH
