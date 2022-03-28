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
#ifndef DATA_ANALYSIS_FUNCTIONAL_2D_HH
#define DATA_ANALYSIS_FUNCTIONAL_2D_HH

#include <cmath>
#include <limits>

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "core/blockStatistics.h"
#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "dataProcessors/dataAnalysisFunctional2D.h"
#include "finiteDifference/fdStencils1D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for BlockLattice ******* */

template <typename T, template <typename U> class Descriptor>
BoxSumRhoBarFunctional2D<T, Descriptor>::BoxSumRhoBarFunctional2D() :
    sumRhoBarId(this->getStatistics().subscribeSum())
{ }

template <typename T, template <typename U> class Descriptor>
void BoxSumRhoBarFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> const &cell = lattice.get(iX, iY);
            statistics.gatherSum(sumRhoBarId, cell.getDynamics().computeRhoBar(cell));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxSumRhoBarFunctional2D<T, Descriptor> *BoxSumRhoBarFunctional2D<T, Descriptor>::clone() const
{
    return new BoxSumRhoBarFunctional2D(*this);
}

template <typename T, template <typename U> class Descriptor>
T BoxSumRhoBarFunctional2D<T, Descriptor>::getSumRhoBar() const
{
    return this->getStatistics().getSum(sumRhoBarId);
}

template <typename T, template <typename U> class Descriptor>
DotSumRhoBarFunctional2D<T, Descriptor>::DotSumRhoBarFunctional2D() :
    sumRhoBarId(this->getStatistics().subscribeSum())
{ }

template <typename T, template <typename U> class Descriptor>
void DotSumRhoBarFunctional2D<T, Descriptor>::process(
    DotList2D const &dotList, BlockLattice2D<T, Descriptor> &lattice)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint index = 0; index < dotList.getN(); ++index) {
        Cell<T, Descriptor> const &cell =
            lattice.get(dotList.getDot(index).x, dotList.getDot(index).y);
        statistics.gatherSum(sumRhoBarId, cell.getDynamics().computeRhoBar(cell));
    }
    this->getStatistics().evaluate();
}

template <typename T, template <typename U> class Descriptor>
DotSumRhoBarFunctional2D<T, Descriptor> *DotSumRhoBarFunctional2D<T, Descriptor>::clone() const
{
    return new DotSumRhoBarFunctional2D(*this);
}

template <typename T, template <typename U> class Descriptor>
T DotSumRhoBarFunctional2D<T, Descriptor>::getSumRhoBar() const
{
    // cout << this->getStatistics().getSum(sumRhoBarId) << endl;
    return this->getStatistics().getSum(sumRhoBarId);
}

template <typename T, template <typename U> class Descriptor>
BoxSumEnergyFunctional2D<T, Descriptor>::BoxSumEnergyFunctional2D() :
    sumEnergyId(this->getStatistics().subscribeSum())
{ }

template <typename T, template <typename U> class Descriptor>
void BoxSumEnergyFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Array<T, Descriptor<T>::d> velocity;
            lattice.get(iX, iY).computeVelocity(velocity);
            T uNormSqr = VectorTemplate<T, Descriptor>::normSqr(velocity);
            statistics.gatherSum(sumEnergyId, uNormSqr);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxSumEnergyFunctional2D<T, Descriptor> *BoxSumEnergyFunctional2D<T, Descriptor>::clone() const
{
    return new BoxSumEnergyFunctional2D(*this);
}

template <typename T, template <typename U> class Descriptor>
T BoxSumEnergyFunctional2D<T, Descriptor>::getSumEnergy() const
{
    return this->getStatistics().getSum(sumEnergyId) / (T)2;
}

template <typename T, template <typename U> class Descriptor>
void BoxSumEnergyFunctional2D<T, Descriptor>::getDimensionsX(std::vector<int> &dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = 2;
}

template <typename T, template <typename U> class Descriptor>
void BoxSumEnergyFunctional2D<T, Descriptor>::getDimensionsT(std::vector<int> &dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = -2;
}

template <typename T, template <typename U> class Descriptor>
DotSumVelocityComponentFunctional2D<T, Descriptor>::DotSumVelocityComponentFunctional2D(
    int iComponent_) :
    sumVelocityComponentId(this->getStatistics().subscribeSum()), iComponent(iComponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void DotSumVelocityComponentFunctional2D<T, Descriptor>::process(
    DotList2D const &dotList, BlockLattice2D<T, Descriptor> &lattice)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint index = 0; index < dotList.getN(); ++index) {
        Cell<T, Descriptor> const &cell =
            lattice.get(dotList.getDot(index).x, dotList.getDot(index).y);
        Array<T, Descriptor<T>::d> velocity(0.0, 0.0);
        cell.getDynamics().computeVelocity(cell, velocity);
        statistics.gatherSum(sumVelocityComponentId, velocity[iComponent]);
    }
    this->getStatistics().evaluate();
}

template <typename T, template <typename U> class Descriptor>
DotSumVelocityComponentFunctional2D<T, Descriptor>
    *DotSumVelocityComponentFunctional2D<T, Descriptor>::clone() const
{
    return new DotSumVelocityComponentFunctional2D(*this);
}

template <typename T, template <typename U> class Descriptor>
T DotSumVelocityComponentFunctional2D<T, Descriptor>::getSumVelocityComponent() const
{
    return this->getStatistics().getSum(sumVelocityComponentId);
}
/* *************** Data Functionals for BlockLattice ***************** */

template <typename T, template <typename U> class Descriptor>
void CopyPopulationsFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &latticeFrom,
    BlockLattice2D<T, Descriptor> &latticeTo)
{
    Dot2D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            latticeTo.get(iX + offset.x, iY + offset.y).attributeValues(latticeFrom.get(iX, iY));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
CopyPopulationsFunctional2D<T, Descriptor> *CopyPopulationsFunctional2D<T, Descriptor>::clone()
    const
{
    return new CopyPopulationsFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void CopyPopulationsFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void CopyConvertPopulationsFunctional2D<T1, Descriptor1, T2, Descriptor2>::process(
    Box2D domain, BlockLattice2D<T1, Descriptor1> &latticeFrom,
    BlockLattice2D<T2, Descriptor2> &latticeTo)
{
    Dot2D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint ipop = 0; ipop < Descriptor1<T1>::q; ++ipop)
                latticeTo.get(iX + offset.x, iY + offset.y)[ipop] =
                    (T2)latticeFrom.get(iX, iY)[ipop];
        }
    }
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
CopyConvertPopulationsFunctional2D<T1, Descriptor1, T2, Descriptor2>
    *CopyConvertPopulationsFunctional2D<T1, Descriptor1, T2, Descriptor2>::clone() const
{
    return new CopyConvertPopulationsFunctional2D<T1, Descriptor1, T2, Descriptor2>(*this);
}

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
void CopyConvertPopulationsFunctional2D<T1, Descriptor1, T2, Descriptor2>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
void LatticeCopyAllFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &latticeFrom,
    BlockLattice2D<T, Descriptor> &latticeTo)
{
    std::vector<char> data;
    Dot2D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> const &fromCell = latticeFrom.get(iX, iY);
            Cell<T, Descriptor> &toCell = latticeTo.get(iX + offset.x, iY + offset.y);
            toCell.attributeValues(fromCell);

            data.clear();
            HierarchicSerializer serializer(data, fromCell.getDynamics().getId());
            fromCell.getDynamics().serialize(serializer);

            HierarchicUnserializer unSerializer(data, 0);
            toCell.getDynamics().unserialize(unSerializer);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeCopyAllFunctional2D<T, Descriptor> *LatticeCopyAllFunctional2D<T, Descriptor>::clone() const
{
    return new LatticeCopyAllFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void LatticeCopyAllFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::allVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT LatticeCopyAllFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void LatticeRegenerateFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &latticeFrom,
    BlockLattice2D<T, Descriptor> &latticeTo)
{
    Dot2D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            latticeTo.get(iX + offset.x, iY + offset.y).attributeValues(latticeFrom.get(iX, iY));
            latticeTo.attributeDynamics(
                iX + offset.x, iY + offset.y, latticeFrom.get(iX, iY).getDynamics().clone());
            latticeTo.get(iX + offset.x, iY + offset.y).attributeValues(latticeFrom.get(iX, iY));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
LatticeRegenerateFunctional2D<T, Descriptor> *LatticeRegenerateFunctional2D<T, Descriptor>::clone()
    const
{
    return new LatticeRegenerateFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void LatticeRegenerateFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    // Full dynamics object must be recreated, because this data processor
    //   re-attributes a new dynamics and acts on the bulk only.
    modified[1] = modif::dataStructure;
}

template <typename T, template <typename U> class Descriptor>
void BoxDensityFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            scalarField.get(iX + offset.x, iY + offset.y) = lattice.get(iX, iY).computeDensity();
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxDensityFunctional2D<T, Descriptor> *BoxDensityFunctional2D<T, Descriptor>::clone() const
{
    return new BoxDensityFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxDensityFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxDensityFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxRhoBarFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> const &cell = lattice.get(iX, iY);
            scalarField.get(iX + offset.x, iY + offset.y) = cell.getDynamics().computeRhoBar(cell);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxRhoBarFunctional2D<T, Descriptor> *BoxRhoBarFunctional2D<T, Descriptor>::clone() const
{
    return new BoxRhoBarFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxRhoBarFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxRhoBarFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxRhoBarJfunctional2D<T, Descriptor>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> fields)
{
    BlockLattice2D<T, Descriptor> &lattice =
        *dynamic_cast<BlockLattice2D<T, Descriptor> *>(fields[0]);
    ScalarField2D<T> &rhoBarField = *dynamic_cast<ScalarField2D<T> *>(fields[1]);
    TensorField2D<T, 2> &jField = *dynamic_cast<TensorField2D<T, 2> *>(fields[2]);
    Dot2D offset1 = computeRelativeDisplacement(lattice, rhoBarField);
    Dot2D offset2 = computeRelativeDisplacement(lattice, jField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> const &cell = lattice.get(iX, iY);
            cell.getDynamics().computeRhoBarJ(
                cell, rhoBarField.get(iX + offset1.x, iY + offset1.y),
                jField.get(iX + offset2.x, iY + offset2.y));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxRhoBarJfunctional2D<T, Descriptor> *BoxRhoBarJfunctional2D<T, Descriptor>::clone() const
{
    return new BoxRhoBarJfunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxRhoBarJfunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // lattice
    modified[1] = modif::staticVariables;  // rhoBar
    modified[2] = modif::staticVariables;  // j
}

template <typename T, template <typename U> class Descriptor>
void PackedRhoBarJfunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, NTensorField2D<T> &rhoBarJField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, rhoBarJField);
    T rhoBar;
    Array<T, 2> j;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> const &cell = lattice.get(iX, iY);
            momentTemplates<T, Descriptor>::get_rhoBar_j(cell, rhoBar, j);
            T *rhoBarJ = rhoBarJField.get(iX + offset.x, iY + offset.y);
            *rhoBarJ = rhoBar;
            j.to_cArray(rhoBarJ + 1);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
PackedRhoBarJfunctional2D<T, Descriptor> *PackedRhoBarJfunctional2D<T, Descriptor>::clone() const
{
    return new PackedRhoBarJfunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void PackedRhoBarJfunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // lattice
    modified[1] = modif::staticVariables;  // rhoBarJ
}

template <typename T, template <typename U> class Descriptor>
void BoxKineticEnergyFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Array<T, Descriptor<T>::d> velocity;
            lattice.get(iX, iY).computeVelocity(velocity);
            scalarField.get(iX + offset.x, iY + offset.y) =
                VectorTemplate<T, Descriptor>::normSqr(velocity) / (T)2;
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxKineticEnergyFunctional2D<T, Descriptor> *BoxKineticEnergyFunctional2D<T, Descriptor>::clone()
    const
{
    return new BoxKineticEnergyFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxKineticEnergyFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxKineticEnergyFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxVelocityNormFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Array<T, Descriptor<T>::d> velocity;
            lattice.get(iX, iY).computeVelocity(velocity);
            // The type cast converts the result of normSqr to type U in case T is of type
            // Complex<U>. Otherwise, the call to std::sqrt would fail, because std::sqrt is
            // overloaded, but not for Palabos' Complex type.
            scalarField.get(iX + offset.x, iY + offset.y) = std::sqrt(
                (typename PlbTraits<T>::BaseType)VectorTemplate<T, Descriptor>::normSqr(velocity));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxVelocityNormFunctional2D<T, Descriptor> *BoxVelocityNormFunctional2D<T, Descriptor>::clone()
    const
{
    return new BoxVelocityNormFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxVelocityNormFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxVelocityNormFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxVelocityComponentFunctional2D<T, Descriptor>::BoxVelocityComponentFunctional2D(int iComponent_) :
    iComponent(iComponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxVelocityComponentFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Array<T, Descriptor<T>::d> velocity;
            lattice.get(iX, iY).computeVelocity(velocity);
            scalarField.get(iX + offset.x, iY + offset.y) = velocity[iComponent];
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxVelocityComponentFunctional2D<T, Descriptor>
    *BoxVelocityComponentFunctional2D<T, Descriptor>::clone() const
{
    return new BoxVelocityComponentFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxVelocityComponentFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxVelocityComponentFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxVelocityFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, Descriptor<T>::d> &tensorField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            lattice.get(iX, iY).computeVelocity(tensorField.get(iX + offset.x, iY + offset.y));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxVelocityFunctional2D<T, Descriptor> *BoxVelocityFunctional2D<T, Descriptor>::clone() const
{
    return new BoxVelocityFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxVelocityFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxVelocityFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxPiNeqFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
{
    Dot2D offset = computeRelativeDisplacement(lattice, PiNeq);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            lattice.get(iX, iY).computePiNeq(PiNeq.get(iX + offset.x, iY + offset.y));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxPiNeqFunctional2D<T, Descriptor> *BoxPiNeqFunctional2D<T, Descriptor>::clone() const
{
    return new BoxPiNeqFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxPiNeqFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxPiNeqFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxShearStressFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq)
{
    Dot2D offset = computeRelativeDisplacement(lattice, PiNeq);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            lattice.get(iX, iY).computeShearStress(PiNeq.get(iX + offset.x, iY + offset.y));
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxShearStressFunctional2D<T, Descriptor> *BoxShearStressFunctional2D<T, Descriptor>::clone() const
{
    return new BoxShearStressFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxShearStressFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxShearStressFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxStrainRateFromStressFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &S)
{
    Dot2D offset = computeRelativeDisplacement(lattice, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> const &cell = lattice.get(iX, iY);
            Array<T, SymmetricTensor<T, Descriptor>::n> &element =
                S.get(iX + offset.x, iY + offset.y);
            cell.computePiNeq(element);
            T omega = cell.getDynamics().getOmega();
            if (cell.getDynamics().hasMoments()) {
                T dynamicOmega =
                    cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                if (!util::isZero(dynamicOmega)) {
                    omega = dynamicOmega;
                }
            }
            T rhoBar = cell.getDynamics().computeRhoBar(cell);
            T prefactor = -omega * Descriptor<T>::invCs2 * Descriptor<T>::invRho(rhoBar) / (T)2;
            for (int iTensor = 0; iTensor < SymmetricTensor<T, Descriptor>::n; ++iTensor) {
                element[iTensor] *= prefactor;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxStrainRateFromStressFunctional2D<T, Descriptor>
    *BoxStrainRateFromStressFunctional2D<T, Descriptor>::clone() const
{
    return new BoxStrainRateFromStressFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxStrainRateFromStressFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxStrainRateFromStressFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxTemperatureFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            scalarField.get(iX + offset.x, iY + offset.y) =
                lattice.get(iX, iY).computeTemperature();
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxTemperatureFunctional2D<T, Descriptor> *BoxTemperatureFunctional2D<T, Descriptor>::clone() const
{
    return new BoxTemperatureFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxTemperatureFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxTemperatureFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxSoundSpeedFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            scalarField.get(iX + offset.x, iY + offset.y) =
                std::sqrt(lattice.get(iX, iY).computeTemperature() * Descriptor<T>::invCs2);
            // the speed of sound for a perfect gas is given by :
            // c = std::sqrt(gamma*p/rho) = std::sqrt(gamma*c_l^2*theta)=std::sqrt(theta/c_l^2)
            // => gamma = 1/cl^4
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxSoundSpeedFunctional2D<T, Descriptor> *BoxSoundSpeedFunctional2D<T, Descriptor>::clone() const
{
    return new BoxSoundSpeedFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxSoundSpeedFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxSoundSpeedFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxPopulationFunctional2D<T, Descriptor>::BoxPopulationFunctional2D(plint iComponent_) :
    iComponent(iComponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxPopulationFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            scalarField.get(iX + offset.x, iY + offset.y) = lattice.get(iX, iY)[iComponent];
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxPopulationFunctional2D<T, Descriptor> *BoxPopulationFunctional2D<T, Descriptor>::clone() const
{
    return new BoxPopulationFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxPopulationFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxPopulationFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxEquilibriumFunctional2D<T, Descriptor>::BoxEquilibriumFunctional2D(plint iComponent_) :
    iComponent(iComponent_)
{ }

template <typename T, template <typename U> class Descriptor>
void BoxEquilibriumFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &equilibrium)
{
    Dot2D offset = computeRelativeDisplacement(lattice, equilibrium);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T rhoBar;
            Array<T, Descriptor<T>::d> j;
            Cell<T, Descriptor> const &cell = lattice.get(iX, iY);
            cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);
            T jSqr = normSqr(j);
            equilibrium.get(iX + offset.x, iY + offset.y) =
                cell.computeEquilibrium(iComponent, rhoBar, j, jSqr);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxEquilibriumFunctional2D<T, Descriptor> *BoxEquilibriumFunctional2D<T, Descriptor>::clone() const
{
    return new BoxEquilibriumFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxEquilibriumFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxEquilibriumFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxAllPopulationsFunctional2D<T, Descriptor>::BoxAllPopulationsFunctional2D()
{ }

template <typename T, template <typename U> class Descriptor>
void BoxAllPopulationsFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, Descriptor<T>::q> &tensorField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                tensorField.get(iX + offset.x, iY + offset.y)[iPop] = lattice.get(iX, iY)[iPop];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxAllPopulationsFunctional2D<T, Descriptor> *BoxAllPopulationsFunctional2D<T, Descriptor>::clone()
    const
{
    return new BoxAllPopulationsFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxAllPopulationsFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxAllPopulationsFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
BoxAllPopulationsToLatticeFunctional2D<T, Descriptor>::BoxAllPopulationsToLatticeFunctional2D()
{ }

template <typename T, template <typename U> class Descriptor>
void BoxAllPopulationsToLatticeFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, Descriptor<T>::q> &tensorField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                lattice.get(iX, iY)[iPop] = tensorField.get(oX, iY + offset.y)[iPop];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxAllPopulationsToLatticeFunctional2D<T, Descriptor>
    *BoxAllPopulationsToLatticeFunctional2D<T, Descriptor>::clone() const
{
    return new BoxAllPopulationsToLatticeFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxAllPopulationsToLatticeFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxAllPopulationsToLatticeFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxOmegaFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Cell<T, Descriptor> const &cell = lattice.get(iX, iY);
            scalarField.get(iX + offset.x, iY + offset.y) = cell.getDynamics().getOmega();
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxOmegaFunctional2D<T, Descriptor> *BoxOmegaFunctional2D<T, Descriptor>::clone() const
{
    return new BoxOmegaFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxOmegaFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxOmegaFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalForceFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, Descriptor<T>::d> &tensorField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T *force = lattice.get(iX, iY).getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
            tensorField.get(iX + offset.x, iY + offset.y).from_cArray(force);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxExternalForceFunctional2D<T, Descriptor> *BoxExternalForceFunctional2D<T, Descriptor>::clone()
    const
{
    return new BoxExternalForceFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalForceFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BoxExternalScalarFunctional2D<T, Descriptor>::BoxExternalScalarFunctional2D(int whichScalar_) :
    whichScalar(whichScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalScalarFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            scalarField.get(oX, oY) = *lattice.get(iX, iY).getExternal(whichScalar);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxExternalScalarFunctional2D<T, Descriptor> *BoxExternalScalarFunctional2D<T, Descriptor>::clone()
    const
{
    return new BoxExternalScalarFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalScalarFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BoxExternalVectorFunctional2D<T, Descriptor>::BoxExternalVectorFunctional2D(int vectorBeginsAt_) :
    vectorBeginsAt(vectorBeginsAt_)
{
    PLB_ASSERT(vectorBeginsAt + Descriptor<T>::d <= Descriptor<T>::ExternalField::numScalars);
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalVectorFunctional2D<T, Descriptor>::process(
    Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
    TensorField2D<T, Descriptor<T>::d> &tensorField)
{
    Dot2D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY = iY + offset.y;
            T *tensor = lattice.get(iX, iY).getExternal(vectorBeginsAt);
            tensorField.get(oX, oY).from_cArray(tensor);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
BoxExternalVectorFunctional2D<T, Descriptor> *BoxExternalVectorFunctional2D<T, Descriptor>::clone()
    const
{
    return new BoxExternalVectorFunctional2D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void BoxExternalVectorFunctional2D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT BoxExternalVectorFunctional2D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for ScalarField ******** */

template <typename T>
BoxScalarSumFunctional2D<T>::BoxScalarSumFunctional2D() :
    sumScalarId(this->getStatistics().subscribeSum())
{ }

template <typename T>
void BoxScalarSumFunctional2D<T>::process(Box2D domain, ScalarField2D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            statistics.gatherSum(sumScalarId, (double)scalarField.get(iX, iY));
        }
    }
}

template <typename T>
BoxScalarSumFunctional2D<T> *BoxScalarSumFunctional2D<T>::clone() const
{
    return new BoxScalarSumFunctional2D<T>(*this);
}

template <typename T>
T BoxScalarSumFunctional2D<T>::getSumScalar() const
{
    double doubleSum = this->getStatistics().getSum(sumScalarId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleSum);
    }
    return (T)doubleSum;
}

template <typename T>
MaskedBoxScalarAverageFunctional2D<T>::MaskedBoxScalarAverageFunctional2D(int flag_) :
    averageScalarId(this->getStatistics().subscribeAverage()), flag(flag_)
{ }

template <typename T>
void MaskedBoxScalarAverageFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &scalarField, ScalarField2D<int> &mask)
{
    Dot2D offset = computeRelativeDisplacement(scalarField, mask);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            if (mask.get(iX + offset.x, iY + offset.y) == flag) {
                statistics.gatherAverage(averageScalarId, (double)scalarField.get(iX, iY));
                statistics.incrementStats();
            }
        }
    }
}

template <typename T>
MaskedBoxScalarAverageFunctional2D<T> *MaskedBoxScalarAverageFunctional2D<T>::clone() const
{
    return new MaskedBoxScalarAverageFunctional2D<T>(*this);
}

template <typename T>
T MaskedBoxScalarAverageFunctional2D<T>::getAverageScalar() const
{
    double doubleAverage = this->getStatistics().getAverage(averageScalarId);
    // The average is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleAverage);
    }
    return (T)doubleAverage;
}

template <typename T>
BoxScalarMinFunctional2D<T>::BoxScalarMinFunctional2D() :
    maxScalarId(this->getStatistics().subscribeMax())
{ }

template <typename T>
void BoxScalarMinFunctional2D<T>::process(Box2D domain, ScalarField2D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            // BlockStatistics computes only maximum, no minimum. Therefore,
            //   the relation min(x) = -max(-x) is used.
            statistics.gatherMax(maxScalarId, -(double)scalarField.get(iX, iY));
        }
    }
}

template <typename T>
BoxScalarMinFunctional2D<T> *BoxScalarMinFunctional2D<T>::clone() const
{
    return new BoxScalarMinFunctional2D<T>(*this);
}

template <typename T>
T BoxScalarMinFunctional2D<T>::getMinScalar() const
{
    // The minus sign accounts for the relation min(x) = -max(-x).
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    double doubleMin = -this->getStatistics().getMax(maxScalarId);
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleMin);
    }
    return (T)doubleMin;
}

template <typename T>
MaskedBoxScalarMinFunctional2D<T>::MaskedBoxScalarMinFunctional2D(int flag_) :
    maxScalarId(this->getStatistics().subscribeMax()), flag(flag_)
{ }

template <typename T>
void MaskedBoxScalarMinFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &scalarField, ScalarField2D<int> &mask)
{
    BlockStatistics &statistics = this->getStatistics();
    Dot2D offset = computeRelativeDisplacement(scalarField, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            if (mask.get(iX + offset.x, iY + offset.y) == flag) {
                // BlockStatistics computes only maximum, no minimum. Therefore,
                //   the relation min(x) = -max(-x) is used.
                statistics.gatherMax(maxScalarId, -(double)scalarField.get(iX, iY));
            }
        }
    }
}

template <typename T>
MaskedBoxScalarMinFunctional2D<T> *MaskedBoxScalarMinFunctional2D<T>::clone() const
{
    return new MaskedBoxScalarMinFunctional2D<T>(*this);
}

template <typename T>
T MaskedBoxScalarMinFunctional2D<T>::getMinScalar() const
{
    // The minus sign accounts for the relation min(x) = -max(-x).
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    double doubleMin = -this->getStatistics().getMax(maxScalarId);
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleMin);
    }
    return (T)doubleMin;
}

template <typename T>
BoxScalarMaxFunctional2D<T>::BoxScalarMaxFunctional2D() :
    maxScalarId(this->getStatistics().subscribeMax())
{ }

template <typename T>
void BoxScalarMaxFunctional2D<T>::process(Box2D domain, ScalarField2D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            statistics.gatherMax(maxScalarId, (double)scalarField.get(iX, iY));
        }
    }
}

template <typename T>
BoxScalarMaxFunctional2D<T> *BoxScalarMaxFunctional2D<T>::clone() const
{
    return new BoxScalarMaxFunctional2D<T>(*this);
}

template <typename T>
T BoxScalarMaxFunctional2D<T>::getMaxScalar() const
{
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    double doubleMax = this->getStatistics().getMax(maxScalarId);
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleMax);
    }
    return (T)doubleMax;
}

template <typename T>
MaskedBoxScalarMaxFunctional2D<T>::MaskedBoxScalarMaxFunctional2D(int flag_) :
    maxScalarId(this->getStatistics().subscribeMax()), flag(flag_)
{ }

template <typename T>
void MaskedBoxScalarMaxFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &scalarField, ScalarField2D<int> &mask)
{
    BlockStatistics &statistics = this->getStatistics();
    Dot2D offset = computeRelativeDisplacement(scalarField, mask);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            if (mask.get(iX + offset.x, iY + offset.y) == flag) {
                statistics.gatherMax(maxScalarId, (double)scalarField.get(iX, iY));
            }
        }
    }
}

template <typename T>
MaskedBoxScalarMaxFunctional2D<T> *MaskedBoxScalarMaxFunctional2D<T>::clone() const
{
    return new MaskedBoxScalarMaxFunctional2D<T>(*this);
}

template <typename T>
T MaskedBoxScalarMaxFunctional2D<T>::getMaxScalar() const
{
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    double doubleMax = this->getStatistics().getMax(maxScalarId);
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleMax);
    }
    return (T)doubleMax;
}

template <typename T>
BoundedBoxScalarSumFunctional2D<T>::BoundedBoxScalarSumFunctional2D() :
    sumScalarId(this->getStatistics().subscribeSum())
{ }

template <typename T>
void BoundedBoxScalarSumFunctional2D<T>::processBulk(Box2D domain, ScalarField2D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            statistics.gatherSum(sumScalarId, (double)scalarField.get(iX, iY));
        }
    }
}

template <typename T>
void BoundedBoxScalarSumFunctional2D<T>::processEdge(
    int direction, int orientation, Box2D domain, ScalarField2D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            // Edge nodes have a weight of 0.5, because only 50% of the
            //   cell centered at the node is inside the computational domain.
            statistics.gatherSum(sumScalarId, (double)scalarField.get(iX, iY) / 2.);
        }
    }
}

template <typename T>
void BoundedBoxScalarSumFunctional2D<T>::processCorner(
    int normalX, int normalY, Box2D domain, ScalarField2D<T> &scalarField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            // Corner nodes have a weight of 0.25, because only 25% of the
            //   cell centered at the node is inside the computational domain.
            statistics.gatherSum(sumScalarId, (double)scalarField.get(iX, iY) / 4.);
        }
    }
}

template <typename T>
BoundedBoxScalarSumFunctional2D<T> *BoundedBoxScalarSumFunctional2D<T>::clone() const
{
    return new BoundedBoxScalarSumFunctional2D<T>(*this);
}

template <typename T>
T BoundedBoxScalarSumFunctional2D<T>::getSumScalar() const
{
    double doubleSum = this->getStatistics().getSum(sumScalarId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T)util::roundToInt(doubleSum);
    }
    return (T)doubleSum;
}

template <typename T1, typename T2>
void CopyConvertScalarFunctional2D<T1, T2>::process(
    Box2D domain, ScalarField2D<T1> &field1, ScalarField2D<T2> &field2)
{
    Dot2D offset = computeRelativeDisplacement(field1, field2);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            field2.get(iX + offset.x, iY + offset.y) = (T2)field1.get(iX, iY);
        }
    }
}

template <typename T1, typename T2>
CopyConvertScalarFunctional2D<T1, T2> *CopyConvertScalarFunctional2D<T1, T2>::clone() const
{
    return new CopyConvertScalarFunctional2D<T1, T2>(*this);
}

template <typename T1, typename T2>
void CopyConvertScalarFunctional2D<T1, T2>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T1, typename T2>
BlockDomain::DomainT CopyConvertScalarFunctional2D<T1, T2>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** Data Functionals for scalar-fields **************** */

template <typename T>
void ExtractScalarSubDomainFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &field1, ScalarField2D<T> &field2)
{
    Dot2D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            field2.get(iX + offset.x, iY + offset.y) = field1.get(iX, iY);
        }
    }
}

template <typename T>
ExtractScalarSubDomainFunctional2D<T> *ExtractScalarSubDomainFunctional2D<T>::clone() const
{
    return new ExtractScalarSubDomainFunctional2D<T>(*this);
}

template <typename T>
void ExtractScalarSubDomainFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT ExtractScalarSubDomainFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** compute absolute value functional 2D ************************************* */

template <typename T>
void ComputeAbsoluteValueFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            B.get(iX + offset.x, iY + offset.y) = std::fabs(A.get(iX, iY));
        }
    }
}

template <typename T>
ComputeAbsoluteValueFunctional2D<T> *ComputeAbsoluteValueFunctional2D<T>::clone() const
{
    return new ComputeAbsoluteValueFunctional2D<T>(*this);
}

template <typename T>
void ComputeAbsoluteValueFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT ComputeAbsoluteValueFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** compute sqrt functional 2D ************************************* */

template <typename T>
void ComputeScalarSqrtFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            B.get(iX + offset.x, iY + offset.y) = std::sqrt(A.get(iX, iY));
        }
    }
}

template <typename T>
ComputeScalarSqrtFunctional2D<T> *ComputeScalarSqrtFunctional2D<T>::clone() const
{
    return new ComputeScalarSqrtFunctional2D<T>(*this);
}

template <typename T>
void ComputeScalarSqrtFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT ComputeScalarSqrtFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** compute pow functional 2D ************************************* */

template <typename T>
ComputeScalarPowFunctional2D<T>::ComputeScalarPowFunctional2D(T power_) : power(power_)
{ }

template <typename T>
void ComputeScalarPowFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            B.get(iX + offset.x, iY + offset.y) = std::pow(A.get(iX, iY), power);
        }
    }
}

template <typename T>
ComputeScalarPowFunctional2D<T> *ComputeScalarPowFunctional2D<T>::clone() const
{
    return new ComputeScalarPowFunctional2D<T>(*this);
}

template <typename T>
void ComputeScalarPowFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT ComputeScalarPowFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** compute log functional 2D ************************************* */

template <typename T>
void ComputeScalarLogFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            PLB_ASSERT(A.get(iX, iY) >= T());
            B.get(iX + offset.x, iY + offset.y) = std::log(A.get(iX, iY));
        }
    }
}

template <typename T>
ComputeScalarLogFunctional2D<T> *ComputeScalarLogFunctional2D<T>::clone() const
{
    return new ComputeScalarLogFunctional2D<T>(*this);
}

template <typename T>
void ComputeScalarLogFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT ComputeScalarLogFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** compute sqrt functional 2D ************************************* */

template <typename T, int nDim>
void ComputeTensorSqrtFunctional2D<T, nDim>::process(
    Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iD = 0; iD < nDim; ++iD) {
                PLB_ASSERT(A.get(iX, iY)[iD] >= T());
                B.get(iX + offset.x, iY + offset.y)[iD] = std::sqrt(A.get(iX, iY)[iD]);
            }
        }
    }
}

template <typename T, int nDim>
ComputeTensorSqrtFunctional2D<T, nDim> *ComputeTensorSqrtFunctional2D<T, nDim>::clone() const
{
    return new ComputeTensorSqrtFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void ComputeTensorSqrtFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT ComputeTensorSqrtFunctional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_lt_alpha_functional2D ************************************* */

template <typename T>
A_lt_alpha_functional2D<T>::A_lt_alpha_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_lt_alpha_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<int> &result)
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offset.x, iY + offset.y) = A.get(iX, iY) < alpha ? 1 : 0;
        }
    }
}

template <typename T>
A_lt_alpha_functional2D<T> *A_lt_alpha_functional2D<T>::clone() const
{
    return new A_lt_alpha_functional2D<T>(*this);
}

template <typename T>
void A_lt_alpha_functional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_lt_alpha_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_gt_alpha_functional2D ************************************* */

template <typename T>
A_gt_alpha_functional2D<T>::A_gt_alpha_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_gt_alpha_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<int> &result)
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offset.x, iY + offset.y) = A.get(iX, iY) > alpha ? 1 : 0;
        }
    }
}

template <typename T>
A_gt_alpha_functional2D<T> *A_gt_alpha_functional2D<T>::clone() const
{
    return new A_gt_alpha_functional2D<T>(*this);
}

template <typename T>
void A_gt_alpha_functional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_gt_alpha_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_plus_alpha_functional2D ************************************* */

template <typename T>
A_plus_alpha_functional2D<T>::A_plus_alpha_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_plus_alpha_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result)
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offset.x, iY + offset.y) = A.get(iX, iY) + alpha;
        }
    }
}

template <typename T>
A_plus_alpha_functional2D<T> *A_plus_alpha_functional2D<T>::clone() const
{
    return new A_plus_alpha_functional2D<T>(*this);
}

template <typename T>
void A_plus_alpha_functional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_plus_alpha_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_minus_alpha_functional2D ************************************** */

template <typename T>
A_minus_alpha_functional2D<T>::A_minus_alpha_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_minus_alpha_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result)
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offset.x, iY + offset.y) = A.get(iX, iY) - alpha;
        }
    }
}

template <typename T>
A_minus_alpha_functional2D<T> *A_minus_alpha_functional2D<T>::clone() const
{
    return new A_minus_alpha_functional2D<T>(*this);
}

template <typename T>
void A_minus_alpha_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_minus_alpha_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Alpha_minus_A_functional2D ************************************* */

template <typename T>
Alpha_minus_A_functional2D<T>::Alpha_minus_A_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void Alpha_minus_A_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result)
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offset.x, iY + offset.y) = alpha - A.get(iX, iY);
        }
    }
}

template <typename T>
Alpha_minus_A_functional2D<T> *Alpha_minus_A_functional2D<T>::clone() const
{
    return new Alpha_minus_A_functional2D<T>(*this);
}

template <typename T>
void Alpha_minus_A_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT Alpha_minus_A_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_times_alpha_functional2D ************************************* */

template <typename T>
A_times_alpha_functional2D<T>::A_times_alpha_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_times_alpha_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result)
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offset.x, iY + offset.y) = A.get(iX, iY) * alpha;
        }
    }
}

template <typename T>
A_times_alpha_functional2D<T> *A_times_alpha_functional2D<T>::clone() const
{
    return new A_times_alpha_functional2D<T>(*this);
}

template <typename T>
void A_times_alpha_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_times_alpha_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_times_alpha_functional2D ************************************* */

template <typename T, int nDim>
Tensor_A_times_alpha_functional2D<T, nDim>::Tensor_A_times_alpha_functional2D(T alpha_) :
    alpha(alpha_)
{ }

template <typename T, int nDim>
void Tensor_A_times_alpha_functional2D<T, nDim>::process(
    Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &result)
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iD = 0; iD < nDim; ++iD) {
                result.get(iX + offset.x, iY + offset.y)[iD] = A.get(iX, iY)[iD] * alpha;
            }
        }
    }
}

template <typename T, int nDim>
Tensor_A_times_alpha_functional2D<T, nDim> *Tensor_A_times_alpha_functional2D<T, nDim>::clone()
    const
{
    return new Tensor_A_times_alpha_functional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_times_alpha_functional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_alpha_functional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_dividedBy_alpha_functional2D ************************************* */

template <typename T>
A_dividedBy_alpha_functional2D<T>::A_dividedBy_alpha_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_dividedBy_alpha_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result)
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offset.x, iY + offset.y) = A.get(iX, iY) / alpha;
        }
    }
}

template <typename T>
A_dividedBy_alpha_functional2D<T> *A_dividedBy_alpha_functional2D<T>::clone() const
{
    return new A_dividedBy_alpha_functional2D<T>(*this);
}

template <typename T>
void A_dividedBy_alpha_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_dividedBy_alpha_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Alpha_dividedBy_A_functional2D ************************************* */

template <typename T>
Alpha_dividedBy_A_functional2D<T>::Alpha_dividedBy_A_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void Alpha_dividedBy_A_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result)
{
    Dot2D offset = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offset.x, iY + offset.y) = alpha / A.get(iX, iY);
        }
    }
}

template <typename T>
Alpha_dividedBy_A_functional2D<T> *Alpha_dividedBy_A_functional2D<T>::clone() const
{
    return new Alpha_dividedBy_A_functional2D<T>(*this);
}

template <typename T>
void Alpha_dividedBy_A_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT Alpha_dividedBy_A_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_plus_alpha_inplace_functional2D ************************************* */

template <typename T>
A_plus_alpha_inplace_functional2D<T>::A_plus_alpha_inplace_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_plus_alpha_inplace_functional2D<T>::process(Box2D domain, ScalarField2D<T> &A)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) += alpha;
        }
    }
}

template <typename T>
A_plus_alpha_inplace_functional2D<T> *A_plus_alpha_inplace_functional2D<T>::clone() const
{
    return new A_plus_alpha_inplace_functional2D<T>(*this);
}

template <typename T>
void A_plus_alpha_inplace_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_plus_alpha_inplace_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_minus_alpha_inplace_functional2D ************************************** */

template <typename T>
A_minus_alpha_inplace_functional2D<T>::A_minus_alpha_inplace_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_minus_alpha_inplace_functional2D<T>::process(Box2D domain, ScalarField2D<T> &A)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) -= alpha;
        }
    }
}

template <typename T>
A_minus_alpha_inplace_functional2D<T> *A_minus_alpha_inplace_functional2D<T>::clone() const
{
    return new A_minus_alpha_inplace_functional2D<T>(*this);
}

template <typename T>
void A_minus_alpha_inplace_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_minus_alpha_inplace_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_times_alpha_inplace_functional2D ************************************* */

template <typename T>
A_times_alpha_inplace_functional2D<T>::A_times_alpha_inplace_functional2D(T alpha_) : alpha(alpha_)
{ }

template <typename T>
void A_times_alpha_inplace_functional2D<T>::process(Box2D domain, ScalarField2D<T> &A)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) *= alpha;
        }
    }
}

template <typename T>
A_times_alpha_inplace_functional2D<T> *A_times_alpha_inplace_functional2D<T>::clone() const
{
    return new A_times_alpha_inplace_functional2D<T>(*this);
}

template <typename T>
void A_times_alpha_inplace_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_times_alpha_inplace_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_dividedBy_alpha_inplace_functional2D ************************************* */

template <typename T>
A_dividedBy_alpha_inplace_functional2D<T>::A_dividedBy_alpha_inplace_functional2D(T alpha_) :
    alpha(alpha_)
{ }

template <typename T>
void A_dividedBy_alpha_inplace_functional2D<T>::process(Box2D domain, ScalarField2D<T> &A)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) /= alpha;
        }
    }
}

template <typename T>
A_dividedBy_alpha_inplace_functional2D<T> *A_dividedBy_alpha_inplace_functional2D<T>::clone() const
{
    return new A_dividedBy_alpha_inplace_functional2D<T>(*this);
}

template <typename T>
void A_dividedBy_alpha_inplace_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_dividedBy_alpha_inplace_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_lt_B_functional2D ****************************************** */

template <typename T>
void A_lt_B_functional2D<T>::processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField2D<T> &A = *dynamic_cast<ScalarField2D<T> *>(fields[0]);
    ScalarField2D<T> &B = *dynamic_cast<ScalarField2D<T> *>(fields[1]);
    ScalarField2D<int> &result = *dynamic_cast<ScalarField2D<int> *>(fields[2]);
    Dot2D offsetB = computeRelativeDisplacement(A, B);
    Dot2D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offsetResult.x, iY + offsetResult.y) =
                A.get(iX, iY) < B.get(iX + offsetB.x, iY + offsetB.y) ? 1 : 0;
        }
    }
}

template <typename T>
A_lt_B_functional2D<T> *A_lt_B_functional2D<T>::clone() const
{
    return new A_lt_B_functional2D<T>(*this);
}

template <typename T>
void A_lt_B_functional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_lt_B_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_gt_B_functional2D ****************************************** */

template <typename T>
void A_gt_B_functional2D<T>::processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField2D<T> &A = *dynamic_cast<ScalarField2D<T> *>(fields[0]);
    ScalarField2D<T> &B = *dynamic_cast<ScalarField2D<T> *>(fields[1]);
    ScalarField2D<int> &result = *dynamic_cast<ScalarField2D<int> *>(fields[2]);
    Dot2D offsetB = computeRelativeDisplacement(A, B);
    Dot2D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offsetResult.x, iY + offsetResult.y) =
                A.get(iX, iY) > B.get(iX + offsetB.x, iY + offsetB.y) ? 1 : 0;
        }
    }
}

template <typename T>
A_gt_B_functional2D<T> *A_gt_B_functional2D<T>::clone() const
{
    return new A_gt_B_functional2D<T>(*this);
}

template <typename T>
void A_gt_B_functional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_gt_B_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_plus_B_functional2D ****************************************** */

template <typename T>
void A_plus_B_functional2D<T>::process(Box2D domain, std::vector<ScalarField2D<T> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField2D<T> &A = *fields[0];
    ScalarField2D<T> &B = *fields[1];
    ScalarField2D<T> &result = *fields[2];
    Dot2D offsetB = computeRelativeDisplacement(A, B);
    Dot2D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offsetResult.x, iY + offsetResult.y) =
                A.get(iX, iY) + B.get(iX + offsetB.x, iY + offsetB.y);
        }
    }
}

template <typename T>
A_plus_B_functional2D<T> *A_plus_B_functional2D<T>::clone() const
{
    return new A_plus_B_functional2D<T>(*this);
}

template <typename T>
void A_plus_B_functional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_plus_B_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_minus_B_functional2D ****************************************** */

template <typename T>
void A_minus_B_functional2D<T>::process(Box2D domain, std::vector<ScalarField2D<T> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField2D<T> &A = *fields[0];
    ScalarField2D<T> &B = *fields[1];
    ScalarField2D<T> &result = *fields[2];
    Dot2D offsetB = computeRelativeDisplacement(A, B);
    Dot2D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offsetResult.x, iY + offsetResult.y) =
                A.get(iX, iY) - B.get(iX + offsetB.x, iY + offsetB.y);
        }
    }
}

template <typename T>
A_minus_B_functional2D<T> *A_minus_B_functional2D<T>::clone() const
{
    return new A_minus_B_functional2D<T>(*this);
}

template <typename T>
void A_minus_B_functional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_minus_B_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_times_B_functional2D ****************************************** */

template <typename T>
void A_times_B_functional2D<T>::process(Box2D domain, std::vector<ScalarField2D<T> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField2D<T> &A = *fields[0];
    ScalarField2D<T> &B = *fields[1];
    ScalarField2D<T> &result = *fields[2];
    Dot2D offsetB = computeRelativeDisplacement(A, B);
    Dot2D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offsetResult.x, iY + offsetResult.y) =
                A.get(iX, iY) * B.get(iX + offsetB.x, iY + offsetB.y);
        }
    }
}

template <typename T>
A_times_B_functional2D<T> *A_times_B_functional2D<T>::clone() const
{
    return new A_times_B_functional2D<T>(*this);
}

template <typename T>
void A_times_B_functional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_times_B_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_dividedBy_B_functional2D ****************************************** */

template <typename T>
void A_dividedBy_B_functional2D<T>::process(Box2D domain, std::vector<ScalarField2D<T> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    ScalarField2D<T> &A = *fields[0];
    ScalarField2D<T> &B = *fields[1];
    ScalarField2D<T> &result = *fields[2];
    Dot2D offsetB = computeRelativeDisplacement(A, B);
    Dot2D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offsetResult.x, iY + offsetResult.y) =
                A.get(iX, iY) / B.get(iX + offsetB.x, iY + offsetB.y);
        }
    }
}

template <typename T>
A_dividedBy_B_functional2D<T> *A_dividedBy_B_functional2D<T>::clone() const
{
    return new A_dividedBy_B_functional2D<T>(*this);
}

template <typename T>
void A_dividedBy_B_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT A_dividedBy_B_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_plus_B_inplace_functional2D ****************************************** */

template <typename T>
void A_plus_B_inplace_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) += B.get(iX + offset.x, iY + offset.y);
        }
    }
}

template <typename T>
A_plus_B_inplace_functional2D<T> *A_plus_B_inplace_functional2D<T>::clone() const
{
    return new A_plus_B_inplace_functional2D<T>(*this);
}

template <typename T>
void A_plus_B_inplace_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT A_plus_B_inplace_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_minus_B_inplace_functional2D ****************************************** */

template <typename T>
void A_minus_B_inplace_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) -= B.get(iX + offset.x, iY + offset.y);
        }
    }
}

template <typename T>
A_minus_B_inplace_functional2D<T> *A_minus_B_inplace_functional2D<T>::clone() const
{
    return new A_minus_B_inplace_functional2D<T>(*this);
}

template <typename T>
void A_minus_B_inplace_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT A_minus_B_inplace_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_times_B_inplace_functional2D ****************************************** */

template <typename T>
void A_times_B_inplace_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) *= B.get(iX + offset.x, iY + offset.y);
        }
    }
}

template <typename T>
A_times_B_inplace_functional2D<T> *A_times_B_inplace_functional2D<T>::clone() const
{
    return new A_times_B_inplace_functional2D<T>(*this);
}

template <typename T>
void A_times_B_inplace_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT A_times_B_inplace_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** A_dividedBy_B_inplace_functional2D ****************************************** */

template <typename T>
void A_dividedBy_B_inplace_functional2D<T>::process(
    Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) /= B.get(iX + offset.x, iY + offset.y);
        }
    }
}

template <typename T>
A_dividedBy_B_inplace_functional2D<T> *A_dividedBy_B_inplace_functional2D<T>::clone() const
{
    return new A_dividedBy_B_inplace_functional2D<T>(*this);
}

template <typename T>
void A_dividedBy_B_inplace_functional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT A_dividedBy_B_inplace_functional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

template <typename T, int nDim>
BoxTensorSumFunctional2D<T, nDim>::BoxTensorSumFunctional2D()
{
    for (plint i = 0; i < nDim; i++) {
        sumTensorId[i] = this->getStatistics().subscribeSum();
    }
}

template <typename T, int nDim>
void BoxTensorSumFunctional2D<T, nDim>::process(Box2D domain, TensorField2D<T, nDim> &tensorField)
{
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint i = 0; i < nDim; i++) {
                statistics.gatherSum(sumTensorId[i], (double)tensorField.get(iX, iY)[i]);
            }
        }
    }
}

template <typename T, int nDim>
BoxTensorSumFunctional2D<T, nDim> *BoxTensorSumFunctional2D<T, nDim>::clone() const
{
    return new BoxTensorSumFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
Array<T, nDim> BoxTensorSumFunctional2D<T, nDim>::getSumTensor() const
{
    Array<T, nDim> sum;
    for (plint i = 0; i < nDim; i++) {
        double doubleSum = this->getStatistics().getSum(sumTensorId[i]);
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            doubleSum = util::roundToInt(doubleSum);
        }
        sum[i] = doubleSum;
    }
    return sum;
}

template <typename T1, typename T2, int nDim>
void CopyConvertTensorFunctional2D<T1, T2, nDim>::process(
    Box2D domain, TensorField2D<T1, nDim> &field1, TensorField2D<T2, nDim> &field2)
{
    Dot2D offset = computeRelativeDisplacement(field1, field2);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (int iDim = 0; iDim < nDim; ++iDim) {
                field2.get(iX + offset.x, iY + offset.y)[iDim] = (T2)field1.get(iX, iY)[iDim];
            }
        }
    }
}

template <typename T1, typename T2, int nDim>
CopyConvertTensorFunctional2D<T1, T2, nDim> *CopyConvertTensorFunctional2D<T1, T2, nDim>::clone()
    const
{
    return new CopyConvertTensorFunctional2D<T1, T2, nDim>(*this);
}

template <typename T1, typename T2, int nDim>
void CopyConvertTensorFunctional2D<T1, T2, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T1, typename T2, int nDim>
BlockDomain::DomainT CopyConvertTensorFunctional2D<T1, T2, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void ExtractTensorSubDomainFunctional2D<T, nDim>::process(
    Box2D domain, TensorField2D<T, nDim> &field1, TensorField2D<T, nDim> &field2)
{
    Dot2D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (int iDim = 0; iDim < nDim; ++iDim) {
                field2.get(iX + offset.x, iY + offset.y)[iDim] = field1.get(iX, iY)[iDim];
            }
        }
    }
}

template <typename T, int nDim>
ExtractTensorSubDomainFunctional2D<T, nDim> *ExtractTensorSubDomainFunctional2D<T, nDim>::clone()
    const
{
    return new ExtractTensorSubDomainFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void ExtractTensorSubDomainFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT ExtractTensorSubDomainFunctional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
ExtractTensorComponentFunctional2D<T, nDim>::ExtractTensorComponentFunctional2D(int iComponent_) :
    iComponent(iComponent_)
{
    PLB_ASSERT(iComponent < nDim);
}

template <typename T, int nDim>
void ExtractTensorComponentFunctional2D<T, nDim>::process(
    Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, nDim> &tensorField)
{
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            scalarField.get(iX, iY) = tensorField.get(iX + offset.x, iY + offset.y)[iComponent];
        }
    }
}

template <typename T, int nDim>
ExtractTensorComponentFunctional2D<T, nDim> *ExtractTensorComponentFunctional2D<T, nDim>::clone()
    const
{
    return new ExtractTensorComponentFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void ExtractTensorComponentFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT ExtractTensorComponentFunctional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void ComputeNormFunctional2D<T, nDim>::process(
    Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, nDim> &tensorField)
{
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            scalarField.get(iX, iY) = std::sqrt(VectorTemplateImpl<T, nDim>::normSqr(
                tensorField.get(iX + offset.x, iY + offset.y)));
        }
    }
}

template <typename T, int nDim>
ComputeNormFunctional2D<T, nDim> *ComputeNormFunctional2D<T, nDim>::clone() const
{
    return new ComputeNormFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void ComputeNormFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT ComputeNormFunctional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void ComputeNormSqrFunctional2D<T, nDim>::process(
    Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, nDim> &tensorField)
{
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T normSqr = T();
            for (plint iD = 0; iD < nDim; ++iD) {
                normSqr += tensorField.get(iX + offset.x, iY + offset.y)[iD]
                           * tensorField.get(iX + offset.x, iY + offset.y)[iD];
            }
            scalarField.get(iX, iY) = normSqr;
        }
    }
}

template <typename T, int nDim>
ComputeNormSqrFunctional2D<T, nDim> *ComputeNormSqrFunctional2D<T, nDim>::clone() const
{
    return new ComputeNormSqrFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void ComputeNormSqrFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT ComputeNormSqrFunctional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void ComputeSymmetricTensorNormFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, 3> &tensorField)
{
    typedef SymmetricTensorImpl<T, 2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Array<T, 3> &el = tensorField.get(iX + offset.x, iY + offset.y);
            scalarField.get(iX, iY) = std::sqrt(
                // Count diagonal components once ...
                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) +
                // .. and off-diagonal component twice, due to symmetry.
                (T)2 * util::sqr(el[tensor::xy]));
        }
    }
}

template <typename T>
ComputeSymmetricTensorNormFunctional2D<T> *ComputeSymmetricTensorNormFunctional2D<T>::clone() const
{
    return new ComputeSymmetricTensorNormFunctional2D<T>(*this);
}

template <typename T>
void ComputeSymmetricTensorNormFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT ComputeSymmetricTensorNormFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void ComputeSymmetricTensorNormSqrFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, 3> &tensorField)
{
    typedef SymmetricTensorImpl<T, 2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Array<T, 3> &el = tensorField.get(iX + offset.x, iY + offset.y);
            scalarField.get(iX, iY) =
                // Count diagonal components once ...
                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) +
                // .. and off-diagonal components twice, due to symmetry.
                (T)2 * util::sqr(el[tensor::xy]);
        }
    }
}

template <typename T>
ComputeSymmetricTensorNormSqrFunctional2D<T> *ComputeSymmetricTensorNormSqrFunctional2D<T>::clone()
    const
{
    return new ComputeSymmetricTensorNormSqrFunctional2D<T>(*this);
}

template <typename T>
void ComputeSymmetricTensorNormSqrFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT ComputeSymmetricTensorNormSqrFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
void ComputeSymmetricTensorTraceFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, 3> &tensorField)
{
    typedef SymmetricTensorImpl<T, 2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            Array<T, 3> &el = tensorField.get(iX + offset.x, iY + offset.y);
            scalarField.get(iX, iY) = el[tensor::xx] + el[tensor::yy];
        }
    }
}

template <typename T>
ComputeSymmetricTensorTraceFunctional2D<T> *ComputeSymmetricTensorTraceFunctional2D<T>::clone()
    const
{
    return new ComputeSymmetricTensorTraceFunctional2D<T>(*this);
}

template <typename T>
void ComputeSymmetricTensorTraceFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT ComputeSymmetricTensorTraceFunctional2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxLocalMaximumPerComponentFunctional2D<T, nDim>::process(
    Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, nDim> &tensorField)
{
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            scalarField.get(iX, iY) = maxElement(tensorField.get(iX + offset.x, iY + offset.y));
        }
    }
}

template <typename T, int nDim>
BoxLocalMaximumPerComponentFunctional2D<T, nDim>
    *BoxLocalMaximumPerComponentFunctional2D<T, nDim>::clone() const
{
    return new BoxLocalMaximumPerComponentFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxLocalMaximumPerComponentFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T>
void BoxGradientFunctional2D<T>::processBulk(
    Box2D domain, ScalarField2D<T> &phi, TensorField2D<T, 2> &gradient)
{
    Dot2D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            gradient.get(iX2, iY2)[0] = fdDataField::bulkXderiv(phi, iX, iY);
            gradient.get(iX2, iY2)[1] = fdDataField::bulkYderiv(phi, iX, iY);
        }
    }
}

template <typename T>
void BoxGradientFunctional2D<T>::processEdge(
    int direction, int orientation, Box2D domain, ScalarField2D<T> &phi,
    TensorField2D<T, 2> &gradient)
{
    Dot2D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            gradient.get(iX2, iY2)[0] =
                fdDataField::edgeXderiv(phi, direction, orientation, iX, iY);
            gradient.get(iX2, iY2)[1] =
                fdDataField::edgeYderiv(phi, direction, orientation, iX, iY);
        }
    }
}

template <typename T>
void BoxGradientFunctional2D<T>::processCorner(
    int normalX, int normalY, Box2D domain, ScalarField2D<T> &phi, TensorField2D<T, 2> &gradient)
{
    Dot2D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            gradient.get(iX2, iY2)[0] = fdDataField::cornerXderiv(phi, normalX, normalY, iX, iY);
            gradient.get(iX2, iY2)[1] = fdDataField::cornerYderiv(phi, normalX, normalY, iX, iY);
        }
    }
}

template <typename T>
BoxGradientFunctional2D<T> *BoxGradientFunctional2D<T>::clone() const
{
    return new BoxGradientFunctional2D<T>(*this);
}

template <typename T>
void BoxGradientFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxGradientFunctional2D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxBulkVorticityFunctional2D<T, nDim>::process(
    Box2D domain, ScalarField2D<T> &vorticity, TensorField2D<T, nDim> &velocity)
{
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            vorticity.get(iX, iY) = fdDataField::bulkVorticity(velocity, iX2, iY2);
        }
    }
}

template <typename T, int nDim>
BoxBulkVorticityFunctional2D<T, nDim> *BoxBulkVorticityFunctional2D<T, nDim>::clone() const
{
    return new BoxBulkVorticityFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxBulkVorticityFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxBulkVorticityFunctional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxVorticityFunctional2D<T, nDim>::processBulk(
    Box2D domain, ScalarField2D<T> &vorticity, TensorField2D<T, nDim> &velocity)
{
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            vorticity.get(iX, iY) = fdDataField::bulkVorticity(velocity, iX2, iY2);
        }
    }
}

template <typename T, int nDim>
void BoxVorticityFunctional2D<T, nDim>::processEdge(
    int direction, int orientation, Box2D domain, ScalarField2D<T> &vorticity,
    TensorField2D<T, nDim> &velocity)
{
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            vorticity.get(iX, iY) =
                fdDataField::edgeVorticity(velocity, direction, orientation, iX2, iY2);
        }
    }
}

template <typename T, int nDim>
void BoxVorticityFunctional2D<T, nDim>::processCorner(
    int normalX, int normalY, Box2D domain, ScalarField2D<T> &vorticity,
    TensorField2D<T, nDim> &velocity)
{
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            vorticity.get(iX, iY) =
                fdDataField::cornerVorticity(velocity, normalX, normalY, iX2, iY2);
        }
    }
}

template <typename T, int nDim>
BoxVorticityFunctional2D<T, nDim> *BoxVorticityFunctional2D<T, nDim>::clone() const
{
    return new BoxVorticityFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxVorticityFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxVorticityFunctional2D<T, nDim>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxBulkStrainRateFunctional2D<T, nDim>::process(
    Box2D domain, TensorField2D<T, nDim> &velocity,
    TensorField2D<T, SymmetricTensorImpl<T, nDim>::n> &S)
{
    typedef SymmetricTensorImpl<T, nDim> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            Array<T, SymmetricTensorImpl<T, nDim>::n> &el = S.get(iX2, iY2);
            el[tensor::xx] = fdDataField::bulkXderiv(velocity, iX, iY, 0);
            el[tensor::yy] = fdDataField::bulkYderiv(velocity, iX, iY, 1);
            el[tensor::xy] = (fdDataField::bulkXderiv(velocity, iX, iY, 1)
                              + fdDataField::bulkYderiv(velocity, iX, iY, 0))
                             / (T)2;
        }
    }
}

template <typename T, int nDim>
BoxBulkStrainRateFunctional2D<T, nDim> *BoxBulkStrainRateFunctional2D<T, nDim>::clone() const
{
    return new BoxBulkStrainRateFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxBulkStrainRateFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxBulkStrainRateFunctional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T, int nDim>
void BoxStrainRateFunctional2D<T, nDim>::processBulk(
    Box2D domain, TensorField2D<T, nDim> &velocity,
    TensorField2D<T, SymmetricTensorImpl<T, nDim>::n> &S)
{
    typedef SymmetricTensorImpl<T, nDim> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            Array<T, SymmetricTensorImpl<T, nDim>::n> &el = S.get(iX2, iY2);
            el[tensor::xx] = fdDataField::bulkXderiv(velocity, iX, iY, 0);
            el[tensor::yy] = fdDataField::bulkYderiv(velocity, iX, iY, 1);
            el[tensor::xy] = (fdDataField::bulkXderiv(velocity, iX, iY, 1)
                              + fdDataField::bulkYderiv(velocity, iX, iY, 0))
                             / (T)2;
        }
    }
}

template <typename T, int nDim>
void BoxStrainRateFunctional2D<T, nDim>::processEdge(
    int direction, int orientation, Box2D domain, TensorField2D<T, nDim> &velocity,
    TensorField2D<T, SymmetricTensorImpl<T, nDim>::n> &S)
{
    typedef SymmetricTensorImpl<T, nDim> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            Array<T, SymmetricTensorImpl<T, nDim>::n> &el = S.get(iX2, iY2);
            el[tensor::xx] = fdDataField::edgeXderiv(velocity, direction, orientation, iX, iY, 0);
            el[tensor::yy] = fdDataField::edgeYderiv(velocity, direction, orientation, iX, iY, 1);
            el[tensor::xy] =
                (fdDataField::edgeXderiv(velocity, direction, orientation, iX, iY, 1)
                 + fdDataField::edgeYderiv(velocity, direction, orientation, iX, iY, 0))
                / (T)2;
        }
    }
}

template <typename T, int nDim>
void BoxStrainRateFunctional2D<T, nDim>::processCorner(
    int normalX, int normalY, Box2D domain, TensorField2D<T, nDim> &velocity,
    TensorField2D<T, SymmetricTensorImpl<T, nDim>::n> &S)
{
    typedef SymmetricTensorImpl<T, nDim> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint iX2 = iX + offset.x;
            plint iY2 = iY + offset.y;
            Array<T, SymmetricTensorImpl<T, nDim>::n> &el = S.get(iX2, iY2);
            el[tensor::xx] = fdDataField::cornerXderiv(velocity, normalX, normalY, iX, iY, 0);
            el[tensor::yy] = fdDataField::cornerYderiv(velocity, normalX, normalY, iX, iY, 1);
            el[tensor::xy] = (fdDataField::cornerXderiv(velocity, normalX, normalY, iX, iY, 1)
                              + fdDataField::cornerYderiv(velocity, normalX, normalY, iX, iY, 0))
                             / (T)2;
        }
    }
}

template <typename T, int nDim>
BoxStrainRateFunctional2D<T, nDim> *BoxStrainRateFunctional2D<T, nDim>::clone() const
{
    return new BoxStrainRateFunctional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void BoxStrainRateFunctional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT BoxStrainRateFunctional2D<T, nDim>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** Tensor_A_plus_B_functional2D ************************************ */

template <typename T, int nDim>
void Tensor_A_plus_B_functional2D<T, nDim>::process(
    Box2D domain, std::vector<TensorField2D<T, nDim> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    TensorField2D<T, nDim> &A = *fields[0];
    TensorField2D<T, nDim> &B = *fields[1];
    TensorField2D<T, nDim> &result = *fields[2];
    Dot2D offsetB = computeRelativeDisplacement(A, B);
    Dot2D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offsetResult.x, iY + offsetResult.y) =
                A.get(iX, iY) + B.get(iX + offsetB.x, iY + offsetB.y);
        }
    }
}

template <typename T, int nDim>
Tensor_A_plus_B_functional2D<T, nDim> *Tensor_A_plus_B_functional2D<T, nDim>::clone() const
{
    return new Tensor_A_plus_B_functional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_plus_B_functional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_plus_B_functional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_minus_B_functional2D ************************************ */

template <typename T, int nDim>
void Tensor_A_minus_B_functional2D<T, nDim>::process(
    Box2D domain, std::vector<TensorField2D<T, nDim> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    TensorField2D<T, nDim> &A = *fields[0];
    TensorField2D<T, nDim> &B = *fields[1];
    TensorField2D<T, nDim> &result = *fields[2];
    Dot2D offsetB = computeRelativeDisplacement(A, B);
    Dot2D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offsetResult.x, iY + offsetResult.y) =
                A.get(iX, iY) - B.get(iX + offsetB.x, iY + offsetB.y);
        }
    }
}

template <typename T, int nDim>
Tensor_A_minus_B_functional2D<T, nDim> *Tensor_A_minus_B_functional2D<T, nDim>::clone() const
{
    return new Tensor_A_minus_B_functional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_minus_B_functional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_minus_B_functional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_times_B_functional2D ************************************ */

template <typename T, int nDim>
void Tensor_A_times_B_functional2D<T, nDim>::process(
    Box2D domain, std::vector<TensorField2D<T, nDim> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    TensorField2D<T, nDim> &A = *fields[0];
    TensorField2D<T, nDim> &B = *fields[1];
    TensorField2D<T, nDim> &result = *fields[2];
    Dot2D offsetB = computeRelativeDisplacement(A, B);
    Dot2D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offsetResult.x, iY + offsetResult.y) =
                A.get(iX, iY) * B.get(iX + offsetB.x, iY + offsetB.y);
        }
    }
}

template <typename T, int nDim>
Tensor_A_times_B_functional2D<T, nDim> *Tensor_A_times_B_functional2D<T, nDim>::clone() const
{
    return new Tensor_A_times_B_functional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_times_B_functional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_B_functional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_dividedBy_B_functional2D ************************************ */

template <typename T, int nDim>
void Tensor_A_dividedBy_B_functional2D<T, nDim>::process(
    Box2D domain, std::vector<TensorField2D<T, nDim> *> fields)
{
    PLB_PRECONDITION(fields.size() == 3);
    TensorField2D<T, nDim> &A = *fields[0];
    TensorField2D<T, nDim> &B = *fields[1];
    TensorField2D<T, nDim> &result = *fields[2];
    Dot2D offsetB = computeRelativeDisplacement(A, B);
    Dot2D offsetResult = computeRelativeDisplacement(A, result);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            result.get(iX + offsetResult.x, iY + offsetResult.y) =
                A.get(iX, iY) / B.get(iX + offsetB.x, iY + offsetB.y);
        }
    }
}

template <typename T, int nDim>
Tensor_A_dividedBy_B_functional2D<T, nDim> *Tensor_A_dividedBy_B_functional2D<T, nDim>::clone()
    const
{
    return new Tensor_A_dividedBy_B_functional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_dividedBy_B_functional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_B_functional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_plus_B_inplace_functional2D ************************************ */

template <typename T, int nDim>
void Tensor_A_plus_B_inplace_functional2D<T, nDim>::process(
    Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) += B.get(iX + offset.x, iY + offset.y);
        }
    }
}

template <typename T, int nDim>
Tensor_A_plus_B_inplace_functional2D<T, nDim>
    *Tensor_A_plus_B_inplace_functional2D<T, nDim>::clone() const
{
    return new Tensor_A_plus_B_inplace_functional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_plus_B_inplace_functional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_plus_B_inplace_functional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_minus_B_inplace_functional2D ************************************ */

template <typename T, int nDim>
void Tensor_A_minus_B_inplace_functional2D<T, nDim>::process(
    Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) -= B.get(iX + offset.x, iY + offset.y);
        }
    }
}

template <typename T, int nDim>
Tensor_A_minus_B_inplace_functional2D<T, nDim>
    *Tensor_A_minus_B_inplace_functional2D<T, nDim>::clone() const
{
    return new Tensor_A_minus_B_inplace_functional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_minus_B_inplace_functional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_minus_B_inplace_functional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_times_B_inplace_functional2D ************************************ */

template <typename T, int nDim>
void Tensor_A_times_B_inplace_functional2D<T, nDim>::process(
    Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) *= B.get(iX + offset.x, iY + offset.y);
        }
    }
}

template <typename T, int nDim>
Tensor_A_times_B_inplace_functional2D<T, nDim>
    *Tensor_A_times_B_inplace_functional2D<T, nDim>::clone() const
{
    return new Tensor_A_times_B_inplace_functional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_times_B_inplace_functional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_B_inplace_functional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** Tensor_A_times_alpha_inplace_functional2D ************************************ */

template <typename T, int nDim>
Tensor_A_times_alpha_inplace_functional2D<T, nDim>::Tensor_A_times_alpha_inplace_functional2D(
    T alpha_) :
    alpha(alpha_)
{ }

template <typename T, int nDim>
void Tensor_A_times_alpha_inplace_functional2D<T, nDim>::process(
    Box2D domain, TensorField2D<T, nDim> &A)
{
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) *= alpha;
        }
    }
}

template <typename T, int nDim>
Tensor_A_times_alpha_inplace_functional2D<T, nDim>
    *Tensor_A_times_alpha_inplace_functional2D<T, nDim>::clone() const
{
    return new Tensor_A_times_alpha_inplace_functional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_times_alpha_inplace_functional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_alpha_inplace_functional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

/* ******** Tensor_A_dividedBy_B_inplace_functional2D ************************************ */

template <typename T, int nDim>
void Tensor_A_dividedBy_B_inplace_functional2D<T, nDim>::process(
    Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B)
{
    Dot2D offset = computeRelativeDisplacement(A, B);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            A.get(iX, iY) /= B.get(iX + offset.x, iY + offset.y);
        }
    }
}

template <typename T, int nDim>
Tensor_A_dividedBy_B_inplace_functional2D<T, nDim>
    *Tensor_A_dividedBy_B_inplace_functional2D<T, nDim>::clone() const
{
    return new Tensor_A_dividedBy_B_inplace_functional2D<T, nDim>(*this);
}

template <typename T, int nDim>
void Tensor_A_dividedBy_B_inplace_functional2D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_B_inplace_functional2D<T, nDim>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONAL_2D_HH
