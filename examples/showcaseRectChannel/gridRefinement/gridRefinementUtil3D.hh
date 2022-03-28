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

/*
 * Several utilities -- implementation.
 */

#ifndef GRID_REFINEMENT_UTIL_3D_HH
#define GRID_REFINEMENT_UTIL_3D_HH

#include <memory>

#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "gridRefinement/couplingInterfaceGenerator3D.h"
#include "gridRefinement/couplingInterfaceGenerator3D.hh"
#include "gridRefinement/gridRefinementUtil3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/multiBlockGenerator3D.hh"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiBlockLattice3D.hh"
#include "multiBlock/nonLocalTransfer3D.h"
#include "multiBlock/nonLocalTransfer3D.hh"
#include "multiGrid/gridConversion3D.h"
#include "multiGrid/gridConversion3D.hh"

namespace plb {

template <
    typename T, template <typename U> class Descriptor,
    template <typename T2, template <typename U2> class Descriptor2> class Engine>
std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > exportLatticeToSingleGridLevel(
    MultiLevelCoupling3D<T, Descriptor, Engine> &lattices, plint exportLevel,
    Box3D const &referenceDomain, plint levelOfReferenceDomain)
{
    PLB_ASSERT(exportLevel >= 0 && exportLevel < lattices.getNumLevels());
    PLB_ASSERT(levelOfReferenceDomain >= 0 && levelOfReferenceDomain < lattices.getNumLevels());

    // First we transform the referenceDomain into coarsest level units to simplify algebra
    // (mainly to avoid divisions completely).
    plint diffCoarsestWithReference = -levelOfReferenceDomain;
    // diffCoarsestWithReference <= 0 always.
    Box3D referenceDomainAtCoarsestLevel =
        referenceDomain.divide(util::intTwoToThePower(-diffCoarsestWithReference));
    Box3D domainAtCoarsestLevel;
#ifdef PLB_DEBUG
    bool intersects =
#endif
        intersect(
            referenceDomainAtCoarsestLevel, lattices.getLevel(0).getBoundingBox(),
            domainAtCoarsestLevel);
    PLB_ASSERT(intersects);

    // Now the coarsest level is the "new reference" level.
    plint diffExportWithCoarsest = exportLevel - 0;
    // diffExportWithCoarsest >= 0 always.
    Box3D domainAtExportLevel =
        domainAtCoarsestLevel.multiply(util::intTwoToThePower(diffExportWithCoarsest));

    /*
    if (exportLevel != 0) {
        if (domainAtExportLevel.x1 == lattices.getLevel(exportLevel).getBoundingBox().x1) {
            domainAtExportLevel.x1--;
        }
        if (domainAtExportLevel.y1 == lattices.getLevel(exportLevel).getBoundingBox().y1) {
            domainAtExportLevel.y1--;
        }
        if (domainAtExportLevel.z1 == lattices.getLevel(exportLevel).getBoundingBox().z1) {
            domainAtExportLevel.z1--;
        }
    }
    */

    std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > singleGridLattice =
        generateMultiBlockLattice(
            domainAtExportLevel, lattices.getLevel(exportLevel).getBackgroundDynamics().clone(),
            lattices.getLevel(exportLevel).getMultiBlockManagement().getEnvelopeWidth());

    for (plint level = 0; level < lattices.getNumLevels(); ++level) {
        // The coarsest level is the "new reference" level.
        plint diffLevelWithCoarsest = level - 0;

        Box3D bbAtLevel = lattices.getLevel(level).getBoundingBox();
        if (level != 0) {
            bbAtLevel.x1--;
            bbAtLevel.y1--;
            bbAtLevel.z1--;
        }
        // diffExportWithCoarsest >= 0 always.
        Box3D domainAtLevel =
            domainAtCoarsestLevel.multiply(util::intTwoToThePower(diffLevelWithCoarsest));

        Box3D bbTmp;
        bool intersects = intersect(domainAtLevel, bbAtLevel, bbTmp);
        if (!intersects) {
            continue;
        }

        // dxScale / dtScale is used to refine / coarsen.
        plint dxScale = level - exportLevel;
        plint dtScale = level - exportLevel;

        std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > subLattice;
        if (dxScale != 0) {
            subLattice = generateMultiBlockLattice(
                bbTmp, lattices.getLevel(level).getBackgroundDynamics().clone(),
                lattices.getLevel(level).getMultiBlockManagement().getEnvelopeWidth());
            copy<T, Descriptor>(
                lattices.getLevel(level), bbTmp, *subLattice, bbTmp, modif::dataStructure);
        }

        if (dxScale > 0) {
            // If dxScale > 0 (<=> level > exportLevel, the export level is lower than the current
            // level) then one needs to coarsen the lattice by a factor 2^(level-exportLevel).

            std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > gridLevel = coarsen(
                *subLattice, dxScale, dtScale,
                lattices.getLevel(level).getBackgroundDynamics().clone());

            // The values of the coarsened lattice are copied in the single grid but only in the
            // closed bound of "level", but in "exportLevel" units.
            Box3D fromDomain(gridLevel->getBoundingBox());
            Box3D toDomain(fromDomain);
            copy<T, Descriptor>(
                *gridLevel, fromDomain, *singleGridLattice, toDomain, modif::dataStructure);
        } else if (dxScale < 0) {
            // If dxScale < 0 (<=>  level < exportLevel, the export level is lower than the current
            // level) then one needs to refine the lattice by a factor 2^(exportLevel-level).
            std::unique_ptr<MultiBlockLattice3D<T, Descriptor> > gridLevel = linearRefine(
                *subLattice, dxScale, dtScale,
                lattices.getLevel(level).getBackgroundDynamics().clone());

            // The values of the refined lattice are copied in the single grid but only in the
            // closed bound of "level", but in "exportLevel" units.
            Box3D fromDomain(gridLevel->getBoundingBox());
            Box3D toDomain(fromDomain);
            copy<T, Descriptor>(
                *gridLevel, fromDomain, *singleGridLattice, toDomain, modif::dataStructure);
        } else {
            // If dxScale == 0 (<=>  level == exportLevel, the export level is the same as the
            // current level) then the values of the lattice are copied in the single grid.
            Box3D fromDomain(bbTmp);
            Box3D toDomain(fromDomain);
            copy<T, Descriptor>(
                lattices.getLevel(level), fromDomain, *singleGridLattice, toDomain,
                modif::dataStructure);
        }
    }

    return singleGridLattice;
}

}  // namespace plb

#endif  // GRID_REFINEMENT_UTIL_3D_HH
