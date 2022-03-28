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

#ifndef MULTI_GRID_DATA_FIELD_3D_HH
#define MULTI_GRID_DATA_FIELD_3D_HH

#include "core/globalDefs.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiGrid/gridConversion3D.h"
#include "multiGrid/multiGridDataField3D.h"
#include "multiGrid/multiGridGenerator3D.h"

namespace plb {

template <typename T>
MultiGridScalarField3D<T>::MultiGridScalarField3D(
    MultiGridManagement3D management_, std::vector<BlockCommunicator3D *> communicators_,
    std::vector<CombinedStatistics *> combinedStatistics_, plint behaviorLevel_) :
    MultiGrid3D(management_, behaviorLevel_)
{
    allocateFields(communicators_, combinedStatistics_);
}

template <typename T>
MultiGridScalarField3D<T>::MultiGridScalarField3D(
    MultiGridManagement3D management_, plint behaviorLevel_) :
    MultiGrid3D(management_, behaviorLevel_)
{
    allocateFields();
}

template <typename T>
MultiGridScalarField3D<T>::MultiGridScalarField3D(MultiGridScalarField3D<T> const &rhs) :
    MultiGrid3D(rhs)
{
    allocateFields();
}

template <typename T>
MultiGridScalarField3D<T>::MultiGridScalarField3D(
    MultiGrid3D const &rhs, Box3D subDomain, bool crop) :
    MultiGrid3D(
        extractManagement(rhs.getMultiGridManagement(), subDomain, crop), rhs.getBehaviorLevel())
{
    allocateFields();
}

template <typename T>
MultiGridScalarField3D<T>::~MultiGridScalarField3D()
{
    for (plint iLevel = 0; iLevel < (plint)fields.size(); ++iLevel) {
        delete fields[iLevel];
    }
}

template <typename T>
void MultiGridScalarField3D<T>::reset()
{
    for (plint iLevel = 0; iLevel < (plint)fields.size(); ++iLevel) {
        fields[iLevel]->reset();
    }
}

template <typename T>
T &MultiGridScalarField3D<T>::get(plint iX, plint iY, plint iZ)
{
    return fields[this->getBehaviorLevel()]->get(iX, iY, iZ);
}

template <typename T>
T const &MultiGridScalarField3D<T>::get(plint iX, plint iY, plint iZ) const
{
    return fields[this->getBehaviorLevel()]->get(iX, iY, iZ);
}

template <typename T>
MultiScalarField3D<T> &MultiGridScalarField3D<T>::getComponent(plint level)
{
    PLB_PRECONDITION(level < (plint)fields.size());
    return *fields[level];
}

template <typename T>
MultiScalarField3D<T> const &MultiGridScalarField3D<T>::getComponent(plint level) const
{
    PLB_PRECONDITION(level < (plint)fields.size());
    return *fields[level];
}

// CAUTION: So far dimDx and dimDt should be 0!
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > MultiGridScalarField3D<T>::convertToCoarsest(
    plint dimDx, plint dimDt)
{
    /**** ADDED BY HELEN MORRISON (2014) ****/
    plint lastLevel = this->getNumLevels() - 1;
    std::unique_ptr<MultiScalarField3D<T> > fieldToCoarsen =
        std::unique_ptr<MultiScalarField3D<T> >(
            new MultiScalarField3D<T>(this->getComponent(lastLevel)));

    copyNonLocal<T>(
        this->getComponent(lastLevel),      // from
        *fieldToCoarsen,                    // to
        fieldToCoarsen->getBoundingBox());  // toDomain

    std::unique_ptr<MultiScalarField3D<T> > tmp;

    for (plint iLevel = lastLevel; iLevel >= 1; --iLevel) {
        tmp = coarsen(*fieldToCoarsen, dimDx, dimDt, 1, 1);

        fieldToCoarsen = generateJoinMultiScalarField<T>(*tmp, this->getComponent(iLevel - 1));

        // copy from reshaped iLevel to result
        copyNonLocal<T>(*tmp, *fieldToCoarsen, fieldToCoarsen->getBoundingBox());
        // copy from iLevel-1 to result
        copyNonLocal<T>(
            this->getComponent(iLevel - 1), *fieldToCoarsen, fieldToCoarsen->getBoundingBox());

        tmp.reset();
    }

    return fieldToCoarsen;

    // TODO modify this to use the methods implemented by Jonas
    //      plint levels = this->getNumLevels();
    //
    //      MultiScalarField3D<T> *copy, *tmp;
    //      copy = joinMultiScalarInCoarsest(
    //          *fields[levels-2],*fields[levels-1], dimDx, dimDt);
    //
    //      tmp=copy;
    //
    //      for (plint iLevel=levels-2; iLevel>0; --iLevel){
    //          copy = joinMultiScalarInCoarsest(
    //              *fields[iLevel-1],*copy, dimDx, dimDt );
    //          delete tmp; // erase the old value of copy
    //          tmp = copy; // keep always a pointer over copy
    //      }
    //
    //      return std::unique_ptr<MultiScalarField3D<T> >(copy);
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > MultiGridScalarField3D<T>::convertToFinest(
    plint dimDx, plint dimDt)
{
    /**** ADDED BY HELEN MORRISON (2014) ****/
    plint lastLevel = this->getNumLevels() - 1;
    std::unique_ptr<MultiScalarField3D<T> > fieldToRefine =
        generateMultiScalarField<T>((this->getComponent(0)).getBoundingBox());

    copyNonLocal<T>(
        this->getComponent(0),             // from
        *fieldToRefine,                    // to
        fieldToRefine->getBoundingBox());  // toDomain

    std::unique_ptr<MultiScalarField3D<T> > tmp;

    for (plint iLevel = 0; iLevel <= lastLevel - 1; ++iLevel) {
        tmp = refine(*fieldToRefine, dimDx, dimDt, -1, -1);

        fieldToRefine =
            generateMultiScalarField<T>((this->getComponent(iLevel + 1)).getBoundingBox());

        // copy from reshaped iLevel to result
        copyNonLocal<T>(*tmp, *fieldToRefine, fieldToRefine->getBoundingBox());
        // copy from iLevel+1 to result
        copyNonLocal<T>(
            this->getComponent(iLevel + 1), *fieldToRefine, fieldToRefine->getBoundingBox());

        tmp.reset();
    }

    return fieldToRefine;

    // TODO modify this to use the methods implemented by Jonas
    //      MultiScalarField3D<T> *copy, *tmp;
    //      copy = joinMultiScalarInFinest(
    //          *fields[0],*fields[1], dimDx, dimDt );
    //
    //      tmp=copy;
    //
    //      for (plint iLevel=2; iLevel<(plint)fields.size(); ++iLevel){
    //          copy = joinMultiScalarInFinest(
    //                      *tmp,*fields[iLevel], dimDx, dimDt );
    //          delete tmp; // erase the old value of copy
    //          tmp = copy; // keep always a pointer over copy
    //      }
    //
    //      return std::unique_ptr<MultiScalarField3D<T> >(copy);
}

/**** ADDED BY HELEN MORRISON (2014) ****/
template <typename T>
std::unique_ptr<MultiScalarField3D<T> > MultiGridScalarField3D<T>::convertToLevel(
    plint level, plint dimDx, plint dimDt)
{
    plint lastLevel = this->getNumLevels() - 1;

    if (level == 0)
        return this->convertToCoarsest();
    if (level == lastLevel)
        return this->convertToFinest();

    else {
        std::unique_ptr<MultiScalarField3D<T> > tmp;

        std::unique_ptr<MultiScalarField3D<T> > fieldToRefine =
            generateMultiScalarField<T>((this->getComponent(0)).getBoundingBox());

        copyNonLocal<T>(
            this->getComponent(0),             // from
            *fieldToRefine,                    // to
            fieldToRefine->getBoundingBox());  // toDomain

        for (plint iLevel = 0; iLevel <= level - 1; ++iLevel) {
            tmp = refine(*fieldToRefine, dimDx, dimDt, -1, -1);

            fieldToRefine = generateJoinMultiScalarField<T>(*tmp, this->getComponent(iLevel + 1));

            // copy from reshaped iLevel to result
            copyNonLocal<T>(*tmp, *fieldToRefine, fieldToRefine->getBoundingBox());

            // copy from iLevel+1 to result
            copyNonLocal<T>(
                this->getComponent(iLevel + 1), *fieldToRefine, fieldToRefine->getBoundingBox());

            tmp.reset();
        }
        std::unique_ptr<MultiScalarField3D<T> > fieldToCoarsen =
            std::unique_ptr<MultiScalarField3D<T> >(
                new MultiScalarField3D<T>(this->getComponent(lastLevel)));

        copyNonLocal<T>(
            this->getComponent(lastLevel),      // from
            *fieldToCoarsen,                    // to
            fieldToCoarsen->getBoundingBox());  // toDomain

        for (plint iLevel = lastLevel; iLevel >= level + 1; --iLevel) {
            tmp = coarsen(*fieldToCoarsen, dimDx, dimDt, 1, 1);

            fieldToCoarsen = generateJoinMultiScalarField<T>(*tmp, this->getComponent(iLevel - 1));

            // copy from reshaped iLevel to result
            copyNonLocal<T>(*tmp, *fieldToCoarsen, fieldToCoarsen->getBoundingBox());
            // copy from iLevel-1 to result
            copyNonLocal<T>(
                this->getComponent(iLevel - 1), *fieldToCoarsen, fieldToCoarsen->getBoundingBox());
        }

        tmp.reset();

        std::unique_ptr<MultiScalarField3D<T> > result =
            generateJoinMultiScalarField<T>(*fieldToCoarsen, *fieldToRefine);

        copyNonLocal<T>(*fieldToRefine, *result, result->getBoundingBox());
        copyNonLocal<T>(this->getComponent(level), *result, result->getBoundingBox());
        copyNonLocal<T>(*fieldToCoarsen, *result, result->getBoundingBox());

        fieldToCoarsen.reset();
        fieldToRefine.reset();

        return result;
    }
}

/// Using the management object, create the corresponding fields for each level
template <typename T>
void MultiGridScalarField3D<T>::allocateFields()
{
    fields = generateScalarFields<T>(
        this->getMultiGridManagement(),
        defaultMultiGridPolicy3D().getBlockCommunicator(this->getNumLevels()),
        defaultMultiGridPolicy3D().getCombinedStatistics(this->getNumLevels()));
}

template <typename T>
void MultiGridScalarField3D<T>::allocateFields(
    std::vector<BlockCommunicator3D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics)
{
    fields =
        generateScalarFields<T>(this->getMultiGridManagement(), communicators, combinedStatistics);
}

template <typename T>
int MultiGridScalarField3D<T>::getBlockId() const
{
    return fields[0]->getStaticId();
}

/* ************** MultiGridTensorField3D ******************* */

template <typename T, int nDim>
MultiGridTensorField3D<T, nDim>::MultiGridTensorField3D(
    MultiGridManagement3D management_, std::vector<BlockCommunicator3D *> communicators_,
    std::vector<CombinedStatistics *> combinedStatistics_, plint behaviorLevel_) :
    MultiGrid3D(management_, behaviorLevel_)
{
    allocateFields(communicators_, combinedStatistics_);
}

template <typename T, int nDim>
MultiGridTensorField3D<T, nDim>::MultiGridTensorField3D(
    MultiGridManagement3D management_, plint behaviorLevel_) :
    MultiGrid3D(management_, behaviorLevel_)
{
    allocateFields(
        defaultMultiGridPolicy3D().getBlockCommunicator(this->getNumLevels()),
        defaultMultiGridPolicy3D().getCombinedStatistics(this->getNumLevels()));
}

// TODO: Unused rhs
template <typename T, int nDim>
MultiGridTensorField3D<T, nDim>::MultiGridTensorField3D(MultiGridTensorField3D<T, nDim> const &rhs)
{
    allocateFields();
}

template <typename T, int nDim>
MultiGridTensorField3D<T, nDim>::MultiGridTensorField3D(MultiGrid3D const &rhs) :
    MultiGrid3D(rhs.getMultiGridManagement(), rhs.getBehaviorLevel())
{
    allocateFields();
}

template <typename T, int nDim>
MultiGridTensorField3D<T, nDim>::MultiGridTensorField3D(
    MultiGrid3D const &rhs, Box3D subDomain, bool crop) :
    MultiGrid3D(
        extractManagement(rhs.getMultiGridManagement(), subDomain, crop), rhs.getBehaviorLevel())
{
    allocateFields();
}

template <typename T, int nDim>
MultiGridTensorField3D<T, nDim>::~MultiGridTensorField3D()
{
    for (pluint iLevel = 0; iLevel < fields.size(); ++iLevel) {
        delete fields[iLevel];
    }
}

template <typename T, int nDim>
void MultiGridTensorField3D<T, nDim>::reset()
{
    for (pluint iLevel = 0; iLevel < fields.size(); ++iLevel) {
        fields[iLevel]->reset();
    }
}

template <typename T, int nDim>
Array<T, nDim> &MultiGridTensorField3D<T, nDim>::get(plint iX, plint iY, plint iZ)
{
    return fields[this->getBehaviorLevel()]->get(iX, iY, iZ);
}

template <typename T, int nDim>
Array<T, nDim> const &MultiGridTensorField3D<T, nDim>::get(plint iX, plint iY, plint iZ) const
{
    return fields[this->getBehaviorLevel()]->get(iX, iY, iZ);
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim> &MultiGridTensorField3D<T, nDim>::getComponent(plint level)
{
    PLB_PRECONDITION(level < (plint)fields.size());
    return *fields[level];
}

template <typename T, int nDim>
MultiTensorField3D<T, nDim> const &MultiGridTensorField3D<T, nDim>::getComponent(plint level) const
{
    PLB_PRECONDITION(level < (plint)fields.size());
    return *fields[level];
}

template <typename T, int nDim>
int MultiGridTensorField3D<T, nDim>::getBlockId() const
{
    return fields[this->getBehaviorLevel()]->getStaticId();
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > MultiGridTensorField3D<T, nDim>::convertToCoarsest(
    plint dimDx, plint dimDt)
{
    /**** ADDED BY HELEN MORRISON (2014) ****/

    plint lastLevel = this->getNumLevels() - 1;
    std::unique_ptr<MultiTensorField3D<T, nDim> > fieldToCoarsen =
        std::unique_ptr<MultiTensorField3D<T, nDim> >(
            new MultiTensorField3D<T, nDim>(this->getComponent(lastLevel)));

    copyNonLocal<T, nDim>(
        this->getComponent(lastLevel),      // from
        *fieldToCoarsen,                    // to
        fieldToCoarsen->getBoundingBox());  // toDomain

    std::unique_ptr<MultiTensorField3D<T, nDim> > tmp;

    for (plint iLevel = lastLevel; iLevel >= 1; --iLevel) {
        tmp = coarsen(*fieldToCoarsen, dimDx, dimDt, 1, 1);

        fieldToCoarsen =
            generateJoinMultiTensorField<T, nDim>(*tmp, this->getComponent(iLevel - 1));

        // copy from reshaped iLevel to result
        copyNonLocal<T, nDim>(*tmp, *fieldToCoarsen, fieldToCoarsen->getBoundingBox());
        // copy from iLevel-1 to result
        copyNonLocal<T, nDim>(
            this->getComponent(iLevel - 1), *fieldToCoarsen, fieldToCoarsen->getBoundingBox());

        tmp.reset();
    }

    return fieldToCoarsen;

    // TODO Complete with the functions that Jonas implemented
    //      plint levels = this->getNumLevels();
    //
    //      MultiTensorField3D<T,nDim> *copy, *tmp;
    //      copy = joinMultiTensorInCoarsest(
    //              *fields[levels-2],*fields[levels-1], dimDx, dimDt);
    //
    //      tmp=copy;
    //
    //      for (plint iLevel=levels-2; iLevel>0; --iLevel){
    //          copy = joinMultiTensorInCoarsest(
    //              *fields[iLevel-1],*copy, dimDx, dimDt );
    //          delete tmp; // erase the old value of copy
    //          tmp = copy; // keep always a pointer over copy
    //      }
    //
    //      return std::unique_ptr<MultiTensorField3D<T,nDim> >(copy);
}

template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > MultiGridTensorField3D<T, nDim>::convertToFinest(
    plint dimDx, plint dimDt)
{
    /**** ADDED BY HELEN MORRISON (2014) ****/
    plint lastLevel = this->getNumLevels() - 1;
    std::unique_ptr<MultiTensorField3D<T, nDim> > fieldToRefine =
        generateMultiTensorField<T, nDim>((this->getComponent(0)).getBoundingBox());

    copyNonLocal<T, nDim>(
        this->getComponent(0),             // from
        *fieldToRefine,                    // to
        fieldToRefine->getBoundingBox());  // toDomain

    std::unique_ptr<MultiTensorField3D<T, nDim> > tmp;

    for (plint iLevel = 0; iLevel <= lastLevel - 1; ++iLevel) {
        tmp = refine(*fieldToRefine, dimDx, dimDt, -1, -1);

        fieldToRefine =
            generateMultiTensorField<T, nDim>((this->getComponent(iLevel + 1)).getBoundingBox());

        // copy from reshaped iLevel to result
        copyNonLocal<T, nDim>(*tmp, *fieldToRefine, fieldToRefine->getBoundingBox());
        // copy from iLevel-1 to result
        copyNonLocal<T, nDim>(
            this->getComponent(iLevel + 1), *fieldToRefine, fieldToRefine->getBoundingBox());

        tmp.reset();
    }

    return fieldToRefine;

    // TODO Complete with the functions that Jonas implemented
    //      MultiTensorField3D<T,nDim> *copy, *tmp;
    //      copy = joinMultiTensorInFinest(
    //          *fields[0],*fields[1], dimDx, dimDt );
    //
    //      tmp=copy;
    //
    //      for (plint iLevel=2; iLevel<(plint)fields.size(); ++iLevel){
    //          copy = joinMultiTensorInFinest(
    //                      *tmp,*fields[iLevel], dimDx, dimDt );
    //          delete tmp; // erase the old value of copy
    //          tmp = copy; // keep always a pointer over copy
    //      }
    //
    //      return std::unique_ptr<MultiTensorField3D<T,nDim> >(copy);
}

/**** ADDED BY HELEN MORRISON (2014) ****/
template <typename T, int nDim>
std::unique_ptr<MultiTensorField3D<T, nDim> > MultiGridTensorField3D<T, nDim>::convertToLevel(
    plint level, plint dimDx, plint dimDt)
{
    plint lastLevel = this->getNumLevels() - 1;

    if (level == 0)
        return this->convertToCoarsest();
    if (level == lastLevel)
        return this->convertToFinest();

    else {
        std::unique_ptr<MultiTensorField3D<T, nDim> > tmp;

        std::unique_ptr<MultiTensorField3D<T, nDim> > fieldToRefine =
            generateMultiTensorField<T, nDim>((this->getComponent(0)).getBoundingBox());

        copyNonLocal<T, nDim>(
            this->getComponent(0),             // from
            *fieldToRefine,                    // to
            fieldToRefine->getBoundingBox());  // toDomain

        for (plint iLevel = 0; iLevel <= level - 1; ++iLevel) {
            tmp = refine(*fieldToRefine, dimDx, dimDt, -1, -1);

            fieldToRefine =
                generateJoinMultiTensorField<T, nDim>(*tmp, this->getComponent(iLevel + 1));

            // copy from reshaped iLevel to result
            copyNonLocal<T, nDim>(*tmp, *fieldToRefine, fieldToRefine->getBoundingBox());

            // copy from iLevel+1 to result
            copyNonLocal<T, nDim>(
                this->getComponent(iLevel + 1), *fieldToRefine, fieldToRefine->getBoundingBox());

            tmp.reset();
        }
        std::unique_ptr<MultiTensorField3D<T, nDim> > fieldToCoarsen =
            std::unique_ptr<MultiTensorField3D<T, nDim> >(
                new MultiTensorField3D<T, nDim>(this->getComponent(lastLevel)));

        copyNonLocal<T>(
            this->getComponent(lastLevel),      // from
            *fieldToCoarsen,                    // to
            fieldToCoarsen->getBoundingBox());  // toDomain

        for (plint iLevel = lastLevel; iLevel >= level + 1; --iLevel) {
            tmp = coarsen(*fieldToCoarsen, dimDx, dimDt, 1, 1);

            fieldToCoarsen =
                generateJoinMultiTensorField<T, nDim>(*tmp, this->getComponent(iLevel - 1));

            // copy from reshaped iLevel to result
            copyNonLocal<T, nDim>(*tmp, *fieldToCoarsen, fieldToCoarsen->getBoundingBox());
            // copy from iLevel-1 to result
            copyNonLocal<T, nDim>(
                this->getComponent(iLevel - 1), *fieldToCoarsen, fieldToCoarsen->getBoundingBox());
        }

        tmp.reset();

        std::unique_ptr<MultiTensorField3D<T, nDim> > result =
            generateJoinMultiTensorField<T, nDim>(*fieldToCoarsen, *fieldToRefine);

        copyNonLocal<T, nDim>(*fieldToRefine, *result, result->getBoundingBox());
        copyNonLocal<T, nDim>(this->getComponent(level), *result, result->getBoundingBox());
        copyNonLocal<T, nDim>(*fieldToCoarsen, *result, result->getBoundingBox());

        fieldToCoarsen.reset();
        fieldToRefine.reset();

        return result;
    }
}

/// Using the management object, create the corresponding fields for each level
template <typename T, int nDim>
void MultiGridTensorField3D<T, nDim>::allocateFields()
{
    fields = generateTensorFields<T, nDim>(
        this->getMultiGridManagement(),
        defaultMultiGridPolicy3D().getBlockCommunicator(this->getNumLevels()),
        defaultMultiGridPolicy3D().getCombinedStatistics(this->getNumLevels()));
}

template <typename T, int nDim>
void MultiGridTensorField3D<T, nDim>::allocateFields(
    std::vector<BlockCommunicator3D *> communicators,
    std::vector<CombinedStatistics *> combinedStatistics)
{
    fields = generateTensorFields<T, nDim>(
        this->getMultiGridManagement(), communicators, combinedStatistics);
}

}  // namespace plb

#endif  // MULTI_GRID_DATA_FIELD_3D_HH
