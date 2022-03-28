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
 * Helper functions for data analysis -- header file.
 */

#ifndef MULTI_LEVEL_WRAPPER_3D_HH
#define MULTI_LEVEL_WRAPPER_3D_HH

#include "dataProcessors/dataAnalysisFunctional3D.h"
#include "gridRefinement/dataAnalysisWrapper3D.h"
#include "gridRefinement/multiLevelWrapper3D.h"

namespace plb {

inline bool doesIntersect(const MultiBlockManagement3D &mbManag, const Box3D &bb)
{
    std::map<plint, Box3D> const &bulks = mbManag.getSparseBlockStructure().getBulks();
    std::map<plint, Box3D>::const_iterator it = bulks.begin();
    for (; it != bulks.end(); ++it) {
        if (doesIntersect(it->second, bb)) {
            return true;
        }
    }
    return false;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
Box3D getConsistentDomainThroughLevels(
    const MultiLevelCoupling3D<T, Descriptor, Engine> &lattices, const Box3D &iniDomain,
    plint levelOfDomain)
{
    if (levelOfDomain == 0)
        return iniDomain;

    plint minLevel = 1;
    for (plint iA = 1; iA <= levelOfDomain; ++iA) {
        if (doesIntersect(lattices.getLevel(iA).getMultiBlockManagement(), iniDomain)) {
            minLevel = iA;
            break;
        }
    }
    return iniDomain.divide(util::intTwoToThePower(levelOfDomain - minLevel));
}

template <class MBclass, class MLclass>
std::vector<MBclass *> extractLevel(std::vector<MLclass *> levels, plint iA)
{
    std::vector<MBclass *> mb;
    for (plint iL = 0; iL < (plint)levels.size(); ++iL) {
        mb.push_back(&levels[iL]->getLevel(iA));
    }
    return mb;
}

template <class MLclass>
plint getNumLevels(std::vector<MLclass *> fields)
{
    plint numLevels = fields[0]->getNumLevels();
    for (plint iA = 1; iA < (plint)fields.size(); ++iA) {
        PLB_ASSERT(
            numLevels == fields[iA]->getNumLevels()
            && "Number of levels must be equal among all fields");
    }
    return numLevels;
}

// ====================== BoxProcessingFunctional3D ====================== //
void applyProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevel3D *> &levels, plint numLevels)
{
    std::vector<std::vector<MultiBlock3D *> > blocksPerLevel(numLevels);
    for (plint iA = 0; iA < numLevels; ++iA) {
        for (plint iB = 0; iB < (plint)levels.size(); ++iB) {
            blocksPerLevel[iA].push_back(&(levels[iB]->getLevel(iA)));
        }
    }

    for (plint iA = 0; iA < numLevels; ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(functional->clone(), currentDomain, blocksPerLevel[iA]);
    }
    delete functional;
}

void integrateProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevel3D *> &levels, plint numLevels, plint levelOfIntegration)
{
    std::vector<std::vector<MultiBlock3D *> > blocksPerLevel(numLevels);
    for (plint iA = 0; iA < numLevels; ++iA) {
        for (plint iB = 0; iB < (plint)levels.size(); ++iB) {
            blocksPerLevel[iA].push_back(&(levels[iB]->getLevel(iA)));
        }
    }

    for (plint iA = 0; iA < numLevels; ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, blocksPerLevel[iA], levelOfIntegration);
    }
    delete functional;
}

void integrateProcessingFunctional(
    BoxProcessingFunctional3D *functional, Box3D domain, plint levelOfDomain, MultiLevel3D &actor,
    std::vector<MultiLevel3D *> levels, plint numLevels, plint levelOfIntegration)
{
    std::vector<std::vector<MultiBlock3D *> > blocksPerLevel(numLevels);
    for (plint iA = 0; iA < numLevels; ++iA) {
        for (plint iB = 0; iB < (plint)levels.size(); ++iB) {
            blocksPerLevel[iA].push_back(&(levels[iB]->getLevel(iA)));
        }
    }

    for (plint iA = 0; iA < numLevels; ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, actor.getLevel(iA), blocksPerLevel[iA],
            levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_L ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices)
{
    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(functional->clone(), currentDomain, lattices.getLevel(iA));
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices, plint levelOfIntegration)
{
    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_LL ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor3> class Engine, typename V,
    template <typename W> class Descriptor2,
    template <typename V2, template <typename W2> class Descriptor4> class Engine2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LL<T, Descriptor, V, Descriptor2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices1,
    MultiLevelCoupling3D<V, Descriptor2, Engine2> &lattices2)
{
    PLB_ASSERT(
        lattices1.getNumLevels() == lattices2.getNumLevels()
        && "MultiLevelCoupling3D must have the save number of levels as the MultiLevelCoupling3D");

    for (plint iA = 0; iA < lattices1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, lattices1.getLevel(iA), lattices2.getLevel(iA));
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor3> class Engine, typename V,
    template <typename W> class Descriptor2,
    template <typename V2, template <typename W2> class Descriptor4> class Engine2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LL<T, Descriptor, V, Descriptor2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices1,
    MultiLevelCoupling3D<V, Descriptor2, Engine2> &lattices2, plint levelOfIntegration)
{
    PLB_ASSERT(
        lattices1.getNumLevels() == lattices2.getNumLevels()
        && "MultiLevelCoupling3D must have the save number of levels as the MultiLevelCoupling3D");

    for (plint iA = 0; iA < lattices1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, lattices1.getLevel(iA), lattices2.getLevel(iA),
            levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_LS ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelScalarField3D<T2> &scalarFields)
{
    PLB_ASSERT(lattices.getNumLevels() == scalarFields.getNumLevels() 
        && "MultiLevelScalarField3D must have the save number of levels as the MultiLevelCoupling3D");
    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        applyProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), scalarFields.getLevel(iA));
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelScalarField3D<T2> &scalarFields, plint levelOfIntegration)
{
    PLB_ASSERT(lattices.getNumLevels() == scalarFields.getNumLevels() && 
        "MultiLevelScalarField3D must have the save number of levels as the MultiLevelCoupling3D");

    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), scalarFields.getLevel(iA),
            levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_LS For Output ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelScalarFieldForOutput3D<T2> &scalarFields)
{
    PLB_ASSERT(lattices.getNumLevels() == scalarFields.getNumLevels() 
        && "MultiLevelScalarFieldForOutput3D must have the save number of levels as the MultiLevelCoupling3D");
    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        applyProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), scalarFields.getLevel(iA));
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelScalarFieldForOutput3D<T2> &scalarFields, plint levelOfIntegration)
{
    PLB_ASSERT(lattices.getNumLevels() == scalarFields.getNumLevels() && 
        "MultiLevelScalarFieldForOutput3D must have the save number of levels as the MultiLevelCoupling3D");

    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), scalarFields.getLevel(iA),
            levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_LT ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorField3D<T2, nDim> &tensorFields)
{
    PLB_ASSERT(lattices.getNumLevels() == tensorFields.getNumLevels() && 
        "MultiLevelTensorField3D must have the save number of levels as the MultiLevelCoupling3D");
    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), tensorFields.getLevel(iA));
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorField3D<T2, nDim> &tensorFields, plint levelOfIntegration)
{
    PLB_ASSERT(lattices.getNumLevels() == tensorFields.getNumLevels() && 
        "MultiLevelTensorField3D must have the save number of levels as the MultiLevelCoupling3D");

    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), tensorFields.getLevel(iA),
            levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_LT For Output ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorFieldForOutput3D<T2, nDim> &tensorFields)
{
    PLB_ASSERT(lattices.getNumLevels() == tensorFields.getNumLevels() && 
        "MultiLevelTensorFieldForOutput3D must have the save number of levels as the MultiLevelCoupling3D");
    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), tensorFields.getLevel(iA));
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelCoupling3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorFieldForOutput3D<T2, nDim> &tensorFields, plint levelOfIntegration)
{
    PLB_ASSERT(lattices.getNumLevels() == tensorFields.getNumLevels() && 
        "MultiLevelTensorFieldForOutput3D must have the save number of levels as the MultiLevelCoupling3D");

    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), tensorFields.getLevel(iA),
            levelOfIntegration);
    }
    delete functional;
}

// ======================================================================= //
// ======================== SCALAR FIELD WRAPPERS ======================== //
// ======================================================================= //

// ====================== BoxProcessingFunctional3D_S ====================== //
template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T> &scalarField)
{
    for (plint iA = 0; iA < scalarField.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(functional->clone(), currentDomain, scalarField.getLevel(iA));
    }
    delete functional;
}

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T> &scalarField, plint levelOfIntegration)
{
    for (plint iA = 0; iA < scalarField.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, scalarField.getLevel(iA), levelOfIntegration);
    }
    delete functional;
}

// ================ BoxProcessingFunctional3D_S For Output ================ //
template <typename T>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T> &scalarField)
{
    for (plint iA = 0; iA < scalarField.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(functional->clone(), currentDomain, scalarField.getLevel(iA));
    }
    delete functional;
}

template <typename T>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_S<T> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T> &scalarField, plint levelOfIntegration)
{
    for (plint iA = 0; iA < scalarField.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, scalarField.getLevel(iA), levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_SS ====================== //
template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T1> &scalarField1, MultiLevelScalarField3D<T2> &scalarField2)
{
    PLB_ASSERT(scalarField1.getNumLevels() == scalarField2.getNumLevels() && 
        "MultiLevelScalarField3D must have the save number of levels as the other MultiLevelScalarField3D");
    for (plint iA = 0; iA < scalarField1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, scalarField1.getLevel(iA),
            scalarField2.getLevel(iA));
    }
    delete functional;
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T1> &scalarField1, MultiLevelScalarField3D<T2> &scalarField2,
    plint levelOfIntegration)
{
    PLB_ASSERT(scalarField1.getNumLevels() == scalarField2.getNumLevels() && 
        "MultiLevelScalarField3D must have the save number of levels as the other MultiLevelScalarField3D");

    for (plint iA = 0; iA < scalarField1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, scalarField1.getLevel(iA),
            scalarField2.getLevel(iA), levelOfIntegration);
    }
    delete functional;
}

// ================ BoxProcessingFunctional3D_SS For Output ================== //
template <typename T1, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T1> &scalarField1,
    MultiLevelScalarFieldForOutput3D<T2> &scalarField2)
{
    PLB_ASSERT(scalarField1.getNumLevels() == scalarField2.getNumLevels() && 
        "MultiLevelScalarFieldForOutput3D must have the save number of levels as the other MultiLevelScalarFieldForOutput3D");
    for (plint iA = 0; iA < scalarField1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, scalarField1.getLevel(iA),
            scalarField2.getLevel(iA));
    }
    delete functional;
}

template <typename T1, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_SS<T1, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T1> &scalarField1,
    MultiLevelScalarFieldForOutput3D<T2> &scalarField2, plint levelOfIntegration)
{
    PLB_ASSERT(scalarField1.getNumLevels() == scalarField2.getNumLevels() && 
        "MultiLevelScalarFieldForOutput3D must have the save number of levels as the other MultiLevelScalarFieldForOutput3D");

    for (plint iA = 0; iA < scalarField1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, scalarField1.getLevel(iA),
            scalarField2.getLevel(iA), levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_ST ====================== //
template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T1> &scalar, MultiLevelTensorField3D<T2, nDim> &tensor)
{
    PLB_ASSERT(scalar.getNumLevels() == tensor.getNumLevels() && 
        "MultiLevelScalarField3D must have the save number of levels as the MultiLevelTensorField3D");

    for (plint iA = 0; iA < scalar.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, scalar.getLevel(iA), tensor.getLevel(iA));
    }
    delete functional;
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarField3D<T1> &scalar, MultiLevelTensorField3D<T2, nDim> &tensor,
    plint levelOfIntegration)
{
    PLB_ASSERT(scalar.getNumLevels() == tensor.getNumLevels() && 
        "MultiLevelScalarField3D must have the save number of levels as the MultiLevelTensorField3D");

    for (plint iA = 0; iA < scalar.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, scalar.getLevel(iA), tensor.getLevel(iA),
            levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_ST For Output ====================== //
template <typename T1, typename T2, int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T1> &scalar,
    MultiLevelTensorFieldForOutput3D<T2, nDim> &tensor)
{
    PLB_ASSERT(scalar.getNumLevels() == tensor.getNumLevels() && 
        "MultiLevelScalarFieldForOutput3D must have the save number of levels as the MultiLevelTensorFieldForOutput3D");

    for (plint iA = 0; iA < scalar.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, scalar.getLevel(iA), tensor.getLevel(iA));
    }
    delete functional;
}

template <typename T1, typename T2, int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_ST<T1, T2, nDim> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelScalarFieldForOutput3D<T1> &scalar,
    MultiLevelTensorFieldForOutput3D<T2, nDim> &tensor, plint levelOfIntegration)
{
    PLB_ASSERT(scalar.getNumLevels() == tensor.getNumLevels() && 
        "MultiLevelScalarFieldForOutput3D must have the save number of levels as the MultiLevelTensorFieldForOutput3D");

    for (plint iA = 0; iA < scalar.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, scalar.getLevel(iA), tensor.getLevel(iA),
            levelOfIntegration);
    }
    delete functional;
}

// ====================== ScalarFieldBoxProcessingFunctional3D ====================== //
template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelScalarField3D<T> *> fields)
{
    plint numLevels = getNumLevels(fields);
    for (plint iA = 0; iA < numLevels; ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain,
            extractLevel<MultiScalarField3D<T>, MultiLevelScalarField3D<T> >(fields, iA));
    }
    delete functional;
}

// ================= ScalarFieldBoxProcessingFunctional3D For Output ================ //
template <typename T>
void applyProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelScalarFieldForOutput3D<T> *> fields)
{
    plint numLevels = getNumLevels(fields);
    for (plint iA = 0; iA < numLevels; ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain,
            extractLevel<MultiScalarField3D<T>, MultiLevelScalarField3D<T> >(fields, iA));
    }
    delete functional;
}

template <typename T>
void integrateProcessingFunctional(
    ScalarFieldBoxProcessingFunctional3D<T> *functional, Box3D domain, plint levelOfDomain,
    std::vector<MultiLevelScalarField3D<T> *> fields, plint levelOfIntegration)
{
    plint numLevels = getNumLevels(fields);
    for (plint iA = 0; iA < numLevels; ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain,
            extractLevel<MultiScalarField3D<T>, MultiLevelScalarField3D<T> >(fields, iA),
            levelOfIntegration);
    }
    delete functional;
}

// ======================================================================= //
// ======================== TENSOR FIELD WRAPPERS ======================== //
// ======================================================================= //

// ====================== BoxProcessingFunctional3D_TT ====================== //
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorField3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorField3D<T2, nDim2> &tensorFields2)
{
    PLB_ASSERT(tensorFields1.getNumLevels() == tensorFields2.getNumLevels() && 
        "MultiLevelTensorField3D must have the save number of levels as the other MultiLevelTensorField3D");
    for (plint iA = 0; iA < tensorFields1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, tensorFields1.getLevel(iA),
            tensorFields2.getLevel(iA));
    }
    delete functional;
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorField3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorField3D<T2, nDim2> &tensorFields2, plint levelOfIntegration)
{
    PLB_ASSERT(tensorFields1.getNumLevels() == tensorFields2.getNumLevels() && 
        "MultiLevelTensorField3D must have the save number of levels as the other MultiLevelTensorField3D");

    for (plint iA = 0; iA < tensorFields1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, tensorFields1.getLevel(iA),
            tensorFields2.getLevel(iA), levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_TT For Output ====================== //
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorFieldForOutput3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorFieldForOutput3D<T2, nDim2> &tensorFields2)
{
    PLB_ASSERT(tensorFields1.getNumLevels() == tensorFields2.getNumLevels() && 
        "MultiLevelTensorFieldForOutput3D must have the save number of levels as the other MultiLevelTensorFieldForOutput3D");
    for (plint iA = 0; iA < tensorFields1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, tensorFields1.getLevel(iA),
            tensorFields2.getLevel(iA));
    }
    delete functional;
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorFieldForOutput3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorFieldForOutput3D<T2, nDim2> &tensorFields2, plint levelOfIntegration)
{
    PLB_ASSERT(tensorFields1.getNumLevels() == tensorFields2.getNumLevels() && 
        "MultiLevelTensorFieldForOutput3D must have the save number of levels as the other MultiLevelTensorFieldForOutput3D");

    for (plint iA = 0; iA < tensorFields1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, tensorFields1.getLevel(iA),
            tensorFields2.getLevel(iA), levelOfIntegration);
    }
    delete functional;
}

// ====================== BoundedBoxProcessingFunctional3D_TT ====================== //
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorField3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorField3D<T2, nDim2> &tensorFields2, plint boundaryWidth)
{
    PLB_ASSERT(tensorFields1.getNumLevels() == tensorFields2.getNumLevels() && 
        "MultiLevelTensorField3D must have the save number of levels as the other MultiLevelTensorField3D");
    for (plint iA = 0; iA < tensorFields1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, tensorFields1.getLevel(iA),
            tensorFields2.getLevel(iA), boundaryWidth);
    }
    delete functional;
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorField3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorField3D<T2, nDim2> &tensorFields2, plint boundaryWidth,
    plint levelOfIntegration)
{
    PLB_ASSERT(tensorFields1.getNumLevels() == tensorFields2.getNumLevels() && 
        "MultiLevelTensorField3D must have the save number of levels as the other MultiLevelTensorField3D");

    for (plint iA = 0; iA < tensorFields1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, tensorFields1.getLevel(iA),
            tensorFields2.getLevel(iA), boundaryWidth, levelOfIntegration);
    }
    delete functional;
}

// ====================== BoundedBoxProcessingFunctional3D_TT For Output ====================== //
template <typename T1, int nDim1, typename T2, int nDim2>
void applyProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorFieldForOutput3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorFieldForOutput3D<T2, nDim2> &tensorFields2, plint boundaryWidth)
{
    PLB_ASSERT(tensorFields1.getNumLevels() == tensorFields2.getNumLevels() && 
        "MultiLevelTensorFieldForOutput3D must have the save number of levels as the other MultiLevelTensorFieldForOutput3D");
    for (plint iA = 0; iA < tensorFields1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, tensorFields1.getLevel(iA),
            tensorFields2.getLevel(iA), boundaryWidth);
    }
    delete functional;
}

template <typename T1, int nDim1, typename T2, int nDim2>
void integrateProcessingFunctional(
    BoundedBoxProcessingFunctional3D_TT<T1, nDim1, T2, nDim2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelTensorFieldForOutput3D<T1, nDim1> &tensorFields1,
    MultiLevelTensorFieldForOutput3D<T2, nDim2> &tensorFields2, plint boundaryWidth,
    plint levelOfIntegration)
{
    PLB_ASSERT(tensorFields1.getNumLevels() == tensorFields2.getNumLevels() && 
        "MultiLevelTensorFieldForOutput3D must have the save number of levels as the other MultiLevelTensorFieldForOutput3D");

    for (plint iA = 0; iA < tensorFields1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, tensorFields1.getLevel(iA),
            tensorFields2.getLevel(iA), boundaryWidth, levelOfIntegration);
    }
    delete functional;
}

// ====================== TensorFieldBoxProcessingFunctional3D ====================== //
template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, const Box3D &domain,
    plint levelOfDomain, std::vector<MultiLevelTensorField3D<T, nDim> *> fields)
{
    plint numLevels = getNumLevels<MultiLevelTensorField3D<T, nDim> >(fields);
    for (plint iA = 0; iA < numLevels; ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain,
            extractLevel<MultiTensorField3D<T, nDim>, MultiLevelTensorField3D<T, nDim> >(
                fields, iA));
    }
    delete functional;
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, const Box3D &domain,
    plint levelOfDomain, std::vector<MultiLevelTensorField3D<T, nDim> *> fields,
    plint levelOfIntegration)
{
    plint numLevels = getNumLevels(fields);
    for (plint iA = 0; iA < numLevels; ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain,
            extractLevel<MultiTensorField3D<T, nDim>, MultiLevelTensorField3D<T, nDim> >(
                fields, iA),
            levelOfIntegration);
    }
    delete functional;
}

// ====================== TensorFieldBoxProcessingFunctional3D for output ====================== //
template <typename T, int nDim>
void applyProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, const Box3D &domain,
    plint levelOfDomain, std::vector<MultiLevelTensorFieldForOutput3D<T, nDim> *> fields)
{
    plint numLevels = getNumLevels<MultiLevelTensorFieldForOutput3D<T, nDim> >(fields);
    for (plint iA = 0; iA < numLevels; ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain,
            extractLevel<MultiTensorField3D<T, nDim>, MultiLevelTensorFieldForOutput3D<T, nDim> >(
                fields, iA));
    }
    delete functional;
}

template <typename T, int nDim>
void integrateProcessingFunctional(
    TensorFieldBoxProcessingFunctional3D<T, nDim> *functional, const Box3D &domain,
    plint levelOfDomain, std::vector<MultiLevelTensorFieldForOutput3D<T, nDim> *> fields,
    plint levelOfIntegration)
{
    plint numLevels = getNumLevels(fields);
    for (plint iA = 0; iA < numLevels; ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain,
            extractLevel<MultiTensorField3D<T, nDim>, MultiLevelTensorFieldForOutput3D<T, nDim> >(
                fields, iA),
            levelOfIntegration);
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
Box3D getConsistentDomainThroughLevels(
    const MultiLevelActions3D<T, Descriptor, Engine> &lattices, const Box3D &iniDomain,
    plint levelOfDomain)
{
    if (levelOfDomain == 0)
        return iniDomain;

    plint minLevel = 1;
    for (plint iA = 1; iA <= levelOfDomain; ++iA) {
        if (doesIntersect(lattices.getLevel(iA).getMultiBlockManagement(), iniDomain)) {
            minLevel = iA;
            break;
        }
    }
    return iniDomain.divide(util::intTwoToThePower(levelOfDomain - minLevel));
}
// ====================== BoxProcessingFunctional3D_L ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelActions3D<T, Descriptor, Engine> &lattices)
{
    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(functional->clone(), currentDomain, lattices.getLevel(iA));
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_L<T, Descriptor> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelActions3D<T, Descriptor, Engine> &lattices, plint levelOfIntegration)
{
    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_LL ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor3> class Engine, typename V,
    template <typename W> class Descriptor2,
    template <typename V2, template <typename W2> class Descriptor4> class Engine2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LL<T, Descriptor, V, Descriptor2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelActions3D<T, Descriptor, Engine> &lattices1,
    MultiLevelActions3D<V, Descriptor2, Engine2> &lattices2)
{
    PLB_ASSERT(
        lattices1.getNumLevels() == lattices2.getNumLevels()
        && "MultiLevelActions3D must have the save number of levels as the MultiLevelActions3D");

    for (plint iA = 0; iA < lattices1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, lattices1.getLevel(iA), lattices2.getLevel(iA));
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor3> class Engine, typename V,
    template <typename W> class Descriptor2,
    template <typename V2, template <typename W2> class Descriptor4> class Engine2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LL<T, Descriptor, V, Descriptor2> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelActions3D<T, Descriptor, Engine> &lattices1,
    MultiLevelActions3D<V, Descriptor2, Engine2> &lattices2, plint levelOfIntegration)
{
    PLB_ASSERT(
        lattices1.getNumLevels() == lattices2.getNumLevels()
        && "MultiLevelActions3D must have the save number of levels as the MultiLevelActions3D");

    for (plint iA = 0; iA < lattices1.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, lattices1.getLevel(iA), lattices2.getLevel(iA),
            levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_LS ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelActions3D<T, Descriptor, Engine> &lattices, MultiLevelScalarField3D<T2> &scalarFields)
{
    PLB_ASSERT(
        lattices.getNumLevels() == scalarFields.getNumLevels()
        && "MultiLevelScalarField3D must have the save number of levels as the "
           "MultiLevelActions3D");
    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        applyProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), scalarFields.getLevel(iA));
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LS<T, Descriptor, T2> *functional, Box3D domain, plint levelOfDomain,
    MultiLevelActions3D<T, Descriptor, Engine> &lattices, MultiLevelScalarField3D<T2> &scalarFields,
    plint levelOfIntegration)
{
    PLB_ASSERT(
        lattices.getNumLevels() == scalarFields.getNumLevels()
        && "MultiLevelScalarField3D must have the save number of levels as the "
           "MultiLevelActions3D");

    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), scalarFields.getLevel(iA),
            levelOfIntegration);
    }
    delete functional;
}

// ====================== BoxProcessingFunctional3D_LT ====================== //
template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void applyProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelActions3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorField3D<T2, nDim> &tensorFields)
{
    PLB_ASSERT(
        lattices.getNumLevels() == tensorFields.getNumLevels()
        && "MultiLevelTensorField3D must have the save number of levels as the "
           "MultiLevelActions3D");
    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));

        applyProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), tensorFields.getLevel(iA));
    }
    delete functional;
}

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine, typename T2,
    int nDim>
void integrateProcessingFunctional(
    BoxProcessingFunctional3D_LT<T, Descriptor, T2, nDim> *functional, Box3D domain,
    plint levelOfDomain, MultiLevelActions3D<T, Descriptor, Engine> &lattices,
    MultiLevelTensorField3D<T2, nDim> &tensorFields, plint levelOfIntegration)
{
    PLB_ASSERT(lattices.getNumLevels() == tensorFields.getNumLevels() && 
        "MultiLevelTensorField3D must have the save number of levels as the MultiLevelActions3MultiLevelActions3D");

    for (plint iA = 0; iA < lattices.getNumLevels(); ++iA) {
        Box3D currentDomain = (iA - levelOfDomain >= 0)
                                  ? domain.multiply(util::intTwoToThePower(iA - levelOfDomain))
                                  : domain.divide(util::intTwoToThePower(levelOfDomain - iA));
        integrateProcessingFunctional(
            functional->clone(), currentDomain, lattices.getLevel(iA), tensorFields.getLevel(iA),
            levelOfIntegration);
    }
    delete functional;
}

}  // namespace plb

#endif  // MULTI_LEVEL_WRAPPER_3D_HH
