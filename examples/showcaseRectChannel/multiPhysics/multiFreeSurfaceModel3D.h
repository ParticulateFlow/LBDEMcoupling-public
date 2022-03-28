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

#ifndef MULTI_FREE_SURFACE_MODEL_3D_H
#define MULTI_FREE_SURFACE_MODEL_3D_H

#include <algorithm>

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiPhysics/freeSurfaceInitializer3D.h"
#include "multiPhysics/freeSurfaceModel3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"

namespace plb {

// Couples the velocities of the two fluids. Fluid 1 sees fluid 2, but fluid 2
// does not see fluid 1.
template <typename T, template <typename U> class Descriptor>
class MultiFreeSurfaceOneWayCoupling3D : public BoxProcessingFunctional3D {
public:
    MultiFreeSurfaceOneWayCoupling3D(T interactionStrength_, T rhoDefault1_, T rhoDefault2_) :
        interactionStrength(interactionStrength_),
        rhoDefault1(rhoDefault1_),
        rhoDefault2(rhoDefault2_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual MultiFreeSurfaceOneWayCoupling3D<T, Descriptor> *clone() const
    {
        return new MultiFreeSurfaceOneWayCoupling3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // 1 Fluid.
        modified[1] = modif::nothing;          // 1 rhoBar.
        modified[2] = modif::staticVariables;  // 1 j.
        modified[3] = modif::nothing;          // 1 Mass.
        modified[4] = modif::nothing;          // 1 Volume fraction.
        modified[5] = modif::nothing;          // 1 Flag-status.
        modified[6] = modif::nothing;          // 1 Normal.
        modified[7] = modif::nothing;          // 1 Interface-lists.
        modified[8] = modif::nothing;          // 1 Curvature.
        modified[9] = modif::nothing;          // 1 Outside density.

        modified[10] = modif::nothing;  // 2 Fluid.
        modified[11] = modif::nothing;  // 2 rhoBar.
        modified[12] = modif::nothing;  // 2 j.
        modified[13] = modif::nothing;  // 2 Mass.
        modified[14] = modif::nothing;  // 2 Volume fraction.
        modified[15] = modif::nothing;  // 2 Flag-status.
        modified[16] = modif::nothing;  // 2 Normal.
        modified[17] = modif::nothing;  // 2 Interface-lists.
        modified[18] = modif::nothing;  // 2 Curvature.
        modified[19] = modif::nothing;  // 2 Outside density.
    }

private:
    T interactionStrength;
    T rhoDefault1, rhoDefault2;
};

template <typename T, template <typename U> class Descriptor>
class MultiFreeSurfaceVelocityContinuityCoupling3D : public BoxProcessingFunctional3D {
public:
    MultiFreeSurfaceVelocityContinuityCoupling3D(T rhoDefault1_, T rhoDefault2_) :
        rhoDefault1(rhoDefault1_), rhoDefault2(rhoDefault2_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual MultiFreeSurfaceVelocityContinuityCoupling3D<T, Descriptor> *clone() const
    {
        return new MultiFreeSurfaceVelocityContinuityCoupling3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // 1 Fluid.
        modified[1] = modif::staticVariables;  // 1 rhoBar.
        modified[2] = modif::staticVariables;  // 1 j.
        modified[3] = modif::nothing;          // 1 Mass.
        modified[4] = modif::staticVariables;  // 1 Volume fraction.
        modified[5] = modif::nothing;          // 1 Flag-status.
        modified[6] = modif::nothing;          // 1 Normal.
        modified[7] = modif::nothing;          // 1 Interface-lists.
        modified[8] = modif::nothing;          // 1 Curvature.
        modified[9] = modif::nothing;          // 1 Outside density.

        modified[10] = modif::nothing;          // 2 Fluid.
        modified[11] = modif::staticVariables;  // 2 rhoBar.
        modified[12] = modif::staticVariables;  // 2 j.
        modified[13] = modif::nothing;          // 2 Mass.
        modified[14] = modif::staticVariables;  // 2 Volume fraction.
        modified[15] = modif::nothing;          // 2 Flag-status.
        modified[16] = modif::nothing;          // 2 Normal.
        modified[17] = modif::nothing;          // 2 Interface-lists.
        modified[18] = modif::nothing;          // 2 Curvature.
        modified[19] = modif::nothing;          // 2 Outside density.
    }

private:
    T rhoDefault1, rhoDefault2;
};

template <typename T, template <typename U> class Descriptor>
class MultiFreeSurfaceRepellingForceCoupling3D : public BoxProcessingFunctional3D {
public:
    MultiFreeSurfaceRepellingForceCoupling3D(
        T interactionStrength_, T rhoDefault1_, T rhoDefault2_) :
        interactionStrength(interactionStrength_),
        rhoDefault1(rhoDefault1_),
        rhoDefault2(rhoDefault2_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual MultiFreeSurfaceRepellingForceCoupling3D<T, Descriptor> *clone() const
    {
        return new MultiFreeSurfaceRepellingForceCoupling3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // 1 Fluid.
        modified[1] = modif::staticVariables;  // 1 rhoBar.
        modified[2] = modif::staticVariables;  // 1 j.
        modified[3] = modif::nothing;          // 1 Mass.
        modified[4] = modif::staticVariables;  // 1 Volume fraction.
        modified[5] = modif::nothing;          // 1 Flag-status.
        modified[6] = modif::nothing;          // 1 Normal.
        modified[7] = modif::nothing;          // 1 Interface-lists.
        modified[8] = modif::nothing;          // 1 Curvature.
        modified[9] = modif::nothing;          // 1 Outside density.

        modified[10] = modif::nothing;          // 2 Fluid.
        modified[11] = modif::staticVariables;  // 2 rhoBar.
        modified[12] = modif::staticVariables;  // 2 j.
        modified[13] = modif::nothing;          // 2 Mass.
        modified[14] = modif::staticVariables;  // 2 Volume fraction.
        modified[15] = modif::nothing;          // 2 Flag-status.
        modified[16] = modif::nothing;          // 2 Normal.
        modified[17] = modif::nothing;          // 2 Interface-lists.
        modified[18] = modif::nothing;          // 2 Curvature.
        modified[19] = modif::nothing;          // 2 Outside density.
    }

private:
    T deltaFunction(T r, T h);

private:
    T interactionStrength;
    T rhoDefault1, rhoDefault2;
};

template <typename T, template <typename U> class Descriptor>
class MultiFreeSurfaceComplexCoupling3D : public BoxProcessingFunctional3D {
public:
    MultiFreeSurfaceComplexCoupling3D(T interactionStrength_, T rhoDefault1_, T rhoDefault2_) :
        interactionStrength(interactionStrength_),
        rhoDefault1(rhoDefault1_),
        rhoDefault2(rhoDefault2_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual MultiFreeSurfaceComplexCoupling3D<T, Descriptor> *clone() const
    {
        return new MultiFreeSurfaceComplexCoupling3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // 1 Fluid.
        modified[1] = modif::staticVariables;  // 1 rhoBar.
        modified[2] = modif::staticVariables;  // 1 j.
        modified[3] = modif::nothing;          // 1 Mass.
        modified[4] = modif::staticVariables;  // 1 Volume fraction.
        modified[5] = modif::nothing;          // 1 Flag-status.
        modified[6] = modif::nothing;          // 1 Normal.
        modified[7] = modif::nothing;          // 1 Interface-lists.
        modified[8] = modif::nothing;          // 1 Curvature.
        modified[9] = modif::nothing;          // 1 Outside density.

        modified[10] = modif::nothing;          // 2 Fluid.
        modified[11] = modif::staticVariables;  // 2 rhoBar.
        modified[12] = modif::staticVariables;  // 2 j.
        modified[13] = modif::nothing;          // 2 Mass.
        modified[14] = modif::staticVariables;  // 2 Volume fraction.
        modified[15] = modif::nothing;          // 2 Flag-status.
        modified[16] = modif::nothing;          // 2 Normal.
        modified[17] = modif::nothing;          // 2 Interface-lists.
        modified[18] = modif::nothing;          // 2 Curvature.
        modified[19] = modif::nothing;          // 2 Outside density.
    }

private:
    T deltaFunction(T r, T h);

private:
    T interactionStrength;
    T rhoDefault1, rhoDefault2;
};

template <typename T, template <typename U> class Descriptor>
struct MultiFreeSurfaceFields3D {
    // static const int envelopeWidth = 2;
    static const int envelopeWidth = 3;
    static const int smallEnvelopeWidth = 1;

    MultiFreeSurfaceFields3D(
        SparseBlockStructure3D const &blockStructure, Dynamics<T, Descriptor> *dynamics1_,
        Dynamics<T, Descriptor> *dynamics2_, T rhoDefault1_, T rhoDefault2_, T surfaceTension1_,
        T surfaceTension2_, T contactAngle1_, T contactAngle2_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force1_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force2_, T interactionStrength_,
        ThreadAttribution *threadAttribution = defaultMultiBlockPolicy3D().getThreadAttribution()) :
        dynamics1(dynamics1_),
        dynamics2(dynamics2_),
        rhoDefault1(rhoDefault1_),
        rhoDefault2(rhoDefault2_),
        surfaceTension1(surfaceTension1_),
        surfaceTension2(surfaceTension2_),
        contactAngle1(contactAngle1_),
        contactAngle2(contactAngle2_),
        force1(force1_),
        force2(force2_),
        lattice1(
            MultiBlockManagement3D(blockStructure, threadAttribution->clone(), smallEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(), dynamics1->clone()),
        lattice2(
            MultiBlockManagement3D(blockStructure, threadAttribution->clone(), smallEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(), dynamics2->clone()),
        helperLists1(lattice1),
        helperLists2(lattice2),
        mass1(lattice1),
        mass2(lattice2),
        flag1(
            MultiBlockManagement3D(blockStructure, threadAttribution->clone(), envelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<int>()),
        flag2(
            MultiBlockManagement3D(blockStructure, threadAttribution->clone(), envelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<int>()),
        volumeFraction1((MultiBlock3D &)flag1),
        volumeFraction2((MultiBlock3D &)flag2),
        curvature1((MultiBlock3D &)flag1),
        curvature2((MultiBlock3D &)flag2),
        outsideDensity1((MultiBlock3D &)flag1),
        outsideDensity2((MultiBlock3D &)flag2),
        rhoBar1(lattice1),
        rhoBar2(lattice2),
        j1(lattice1),
        j2(lattice2),
        normal1((MultiBlock3D &)flag1),
        normal2((MultiBlock3D &)flag2),
        interactionStrength(interactionStrength_)
    {
        // MultiFreeSurfaceFields3D does not work with incompressible dynamics at the moment.
        PLB_ASSERT(!dynamics1->velIsJ());
        PLB_ASSERT(!dynamics2->velIsJ());

        delete threadAttribution;

        useSurfaceTension1 = !util::isZero(surfaceTension1);
        useSurfaceTension2 = !util::isZero(surfaceTension2);

        freeSurfaceArgs1 = aggregateFreeSurfaceParams(
            lattice1, rhoBar1, j1, mass1, volumeFraction1, flag1, normal1, helperLists1, curvature1,
            outsideDensity1);

        freeSurfaceArgs2 = aggregateFreeSurfaceParams(
            lattice2, rhoBar2, j2, mass2, volumeFraction2, flag2, normal2, helperLists2, curvature2,
            outsideDensity2);

        multiFreeSurfaceArgs = freeSurfaceArgs1;
        multiFreeSurfaceArgs.insert(
            multiFreeSurfaceArgs.end(), freeSurfaceArgs2.begin(), freeSurfaceArgs2.end());

        initializeInterfaceLists3D<T, Descriptor>(helperLists1);
        initializeInterfaceLists3D<T, Descriptor>(helperLists2);
        lattice1.periodicity().toggleAll(true);
        lattice2.periodicity().toggleAll(true);
        mass1.periodicity().toggleAll(true);
        mass2.periodicity().toggleAll(true);
        flag1.periodicity().toggleAll(true);
        flag2.periodicity().toggleAll(true);
        volumeFraction1.periodicity().toggleAll(true);
        volumeFraction2.periodicity().toggleAll(true);
        curvature1.periodicity().toggleAll(true);
        curvature2.periodicity().toggleAll(true);
        outsideDensity1.periodicity().toggleAll(true);
        outsideDensity2.periodicity().toggleAll(true);

        rhoBar1.periodicity().toggleAll(true);
        rhoBar2.periodicity().toggleAll(true);
        j1.periodicity().toggleAll(true);
        j2.periodicity().toggleAll(true);
        normal1.periodicity().toggleAll(true);
        normal2.periodicity().toggleAll(true);
        // setToConstant(flag1, flag1.getBoundingBox(), (int) freeSurfaceFlag::empty);
        // setToConstant(flag2, flag2.getBoundingBox(), (int) freeSurfaceFlag::empty);
        // setToConstant(outsideDensity1, outsideDensity1.getBoundingBox(), rhoDefault1);
        // setToConstant(outsideDensity2, outsideDensity2.getBoundingBox(), rhoDefault2);
        rhoBarJparam1.push_back(&lattice1);
        rhoBarJparam1.push_back(&rhoBar1);
        rhoBarJparam1.push_back(&j1);
        rhoBarJparam2.push_back(&lattice2);
        rhoBarJparam2.push_back(&rhoBar2);
        rhoBarJparam2.push_back(&j2);

        lattice1.internalStatSubscription().subscribeSum();     // Total mass.
        lattice1.internalStatSubscription().subscribeSum();     // Lost mass.
        lattice1.internalStatSubscription().subscribeIntSum();  // Num interface cells.

        lattice2.internalStatSubscription().subscribeSum();     // Total mass.
        lattice2.internalStatSubscription().subscribeSum();     // Lost mass.
        lattice2.internalStatSubscription().subscribeIntSum();  // Num interface cells.

        freeSurfaceDataProcessors();
    }

    MultiFreeSurfaceFields3D(MultiFreeSurfaceFields3D<T, Descriptor> const &rhs) :
        dynamics1(rhs.dynamics1->clone()),
        dynamics2(rhs.dynamics2->clone()),
        rhoDefault1(rhs.rhoDefault1),
        rhoDefault2(rhs.rhoDefault2),
        surfaceTension1(rhs.surfaceTension1),
        surfaceTension2(rhs.surfaceTension2),
        contactAngle1(rhs.contactAngle1),
        contactAngle2(rhs.contactAngle2),
        useSurfaceTension1(rhs.useSurfaceTension1),
        useSurfaceTension2(rhs.useSurfaceTension2),
        force1(rhs.force1),
        force2(rhs.force2),
        lattice1(rhs.lattice1),
        lattice2(rhs.lattice2),
        helperLists1(rhs.helperLists1),
        helperLists2(rhs.helperLists2),
        mass1(rhs.mass1),
        mass2(rhs.mass2),
        flag1(rhs.flag1),
        flag2(rhs.flag2),
        volumeFraction1(rhs.volumeFraction1),
        volumeFraction2(rhs.volumeFraction2),
        curvature1(rhs.curvature1),
        curvature2(rhs.curvature2),
        outsideDensity1(rhs.outsideDensity1),
        outsideDensity2(rhs.outsideDensity2),
        rhoBar1(rhs.rhoBar1),
        rhoBar2(rhs.rhoBar2),
        j1(rhs.j1),
        j2(rhs.j2),
        normal1(rhs.normal1),
        normal2(rhs.normal2),
        rhoBarJparam1(rhs.rhoBarJparam1),
        rhoBarJparam2(rhs.rhoBarJparam2),
        freeSurfaceArgs1(rhs.freeSurfaceArgs1),
        freeSurfaceArgs2(rhs.freeSurfaceArgs2),
        multiFreeSurfaceArgs(rhs.multiFreeSurfaceArgs),
        interactionStrength(rhs.interactionStrength)
    { }

    void swap(MultiFreeSurfaceFields3D<T, Descriptor> &rhs)
    {
        std::swap(dynamics1, rhs.dynamics1);
        std::swap(dynamics2, rhs.dynamics2);
        std::swap(rhoDefault1, rhs.rhoDefault1);
        std::swap(rhoDefault2, rhs.rhoDefault2);
        std::swap(surfaceTension1, rhs.surfaceTension1);
        std::swap(surfaceTension2, rhs.surfaceTension2);
        std::swap(contactAngle1, rhs.contactAngle1);
        std::swap(contactAngle2, rhs.contactAngle2);
        std::swap(useSurfaceTension1, rhs.useSurfaceTension1);
        std::swap(useSurfaceTension2, rhs.useSurfaceTension2);
        std::swap(force1, rhs.force1);
        std::swap(force2, rhs.force2);
        std::swap(lattice1, rhs.lattice1);
        std::swap(lattice2, rhs.lattice2);
        std::swap(helperLists1, rhs.helperLists1);
        std::swap(helperLists2, rhs.helperLists2);
        std::swap(mass1, rhs.mass1);
        std::swap(mass2, rhs.mass2);
        std::swap(flag1, rhs.flag1);
        std::swap(flag2, rhs.flag2);
        std::swap(volumeFraction1, rhs.volumeFraction1);
        std::swap(volumeFraction2, rhs.volumeFraction2);
        std::swap(curvature1, rhs.curvature1);
        std::swap(curvature2, rhs.curvature2);
        std::swap(outsideDensity1, rhs.outsideDensity1);
        std::swap(outsideDensity2, rhs.outsideDensity2);
        std::swap(rhoBar1, rhs.rhoBar1);
        std::swap(rhoBar2, rhs.rhoBar2);
        std::swap(j1, rhs.j1);
        std::swap(j2, rhs.j2);
        std::swap(normal1, rhs.normal1);
        std::swap(normal2, rhs.normal2);
        std::swap(rhoBarJparam1, rhs.rhoBarJparam1);
        std::swap(rhoBarJparam2, rhs.rhoBarJparam2);
        std::swap(freeSurfaceArgs1, rhs.freeSurfaceArgs1);
        std::swap(freeSurfaceArgs2, rhs.freeSurfaceArgs2);
        std::swap(multiFreeSurfaceArgs, rhs.multiFreeSurfaceArgs);
        std::swap(interactionStrength, rhs.interactionStrength);
    }

    MultiFreeSurfaceFields3D<T, Descriptor> &operator=(
        MultiFreeSurfaceFields3D<T, Descriptor> const &rhs)
    {
        MultiFreeSurfaceFields3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }

    MultiFreeSurfaceFields3D<T, Descriptor> *clone() const
    {
        return new MultiFreeSurfaceFields3D<T, Descriptor>(*this);
    }

    ~MultiFreeSurfaceFields3D()
    {
        delete dynamics1;
        delete dynamics2;
    }

    void periodicityToggle(plint direction, bool periodic)
    {
        PLB_ASSERT(direction == 0 || direction == 1 || direction == 2);

        lattice1.periodicity().toggle(direction, periodic);
        lattice2.periodicity().toggle(direction, periodic);
        mass1.periodicity().toggle(direction, periodic);
        mass2.periodicity().toggle(direction, periodic);
        flag1.periodicity().toggle(direction, periodic);
        flag2.periodicity().toggle(direction, periodic);
        volumeFraction1.periodicity().toggle(direction, periodic);
        volumeFraction2.periodicity().toggle(direction, periodic);
        curvature1.periodicity().toggle(direction, periodic);
        curvature2.periodicity().toggle(direction, periodic);
        outsideDensity1.periodicity().toggle(direction, periodic);
        outsideDensity2.periodicity().toggle(direction, periodic);
        rhoBar1.periodicity().toggle(direction, periodic);
        rhoBar2.periodicity().toggle(direction, periodic);
        j1.periodicity().toggle(direction, periodic);
        j2.periodicity().toggle(direction, periodic);
        normal1.periodicity().toggle(direction, periodic);
        normal2.periodicity().toggle(direction, periodic);
    }

    void periodicityToggleAll(bool periodic)
    {
        lattice1.periodicity().toggleAll(periodic);
        lattice2.periodicity().toggleAll(periodic);
        mass1.periodicity().toggleAll(periodic);
        mass2.periodicity().toggleAll(periodic);
        flag1.periodicity().toggleAll(periodic);
        flag2.periodicity().toggleAll(periodic);
        volumeFraction1.periodicity().toggleAll(periodic);
        volumeFraction2.periodicity().toggleAll(periodic);
        curvature1.periodicity().toggleAll(periodic);
        curvature2.periodicity().toggleAll(periodic);
        outsideDensity1.periodicity().toggleAll(periodic);
        outsideDensity2.periodicity().toggleAll(periodic);
        rhoBar1.periodicity().toggleAll(periodic);
        rhoBar2.periodicity().toggleAll(periodic);
        j1.periodicity().toggleAll(periodic);
        j2.periodicity().toggleAll(periodic);
        normal1.periodicity().toggleAll(periodic);
        normal2.periodicity().toggleAll(periodic);
    }

    void defaultInitialize()
    {
        applyProcessingFunctional(
            new DefaultInitializeFreeSurface3D<T, Descriptor>(
                dynamics1->clone(), force1, rhoDefault1),
            lattice1.getBoundingBox(), freeSurfaceArgs1);

        applyProcessingFunctional(
            new DefaultInitializeFreeSurface3D<T, Descriptor>(
                dynamics2->clone(), force2, rhoDefault2),
            lattice2.getBoundingBox(), freeSurfaceArgs2);
    }

    void partiallyDefaultInitialize()
    {
        applyProcessingFunctional(
            new PartiallyDefaultInitializeFreeSurface3D<T, Descriptor>(
                dynamics1->clone(), force1, rhoDefault1),
            lattice1.getBoundingBox(), freeSurfaceArgs1);

        applyProcessingFunctional(
            new PartiallyDefaultInitializeFreeSurface3D<T, Descriptor>(
                dynamics2->clone(), force2, rhoDefault2),
            lattice2.getBoundingBox(), freeSurfaceArgs2);
    }

    void freeSurfaceDataProcessors()
    {
        MultiBlock3D &actor1 = *freeSurfaceArgs2[0];
        MultiBlock3D &actor2 = *freeSurfaceArgs2[0];

        plint pl;  // Processor level.

        // MultiFreeSurfaceFields3D does not work with incompressible dynamics at the moment.
        bool incompressibleModel1 = false;
        bool incompressibleModel2 = false;

        /***** Initial level ******/
        pl = 0;

        integrateProcessingFunctional(
            new ExternalRhoJcollideAndStream3D<T, Descriptor>, lattice1.getBoundingBox(),
            rhoBarJparam1, pl);
        integrateProcessingFunctional(
            new ExternalRhoJcollideAndStream3D<T, Descriptor>, lattice2.getBoundingBox(),
            rhoBarJparam2, pl);

        integrateProcessingFunctional(
            new FreeSurfaceComputeNormals3D<T, Descriptor>, lattice1.getBoundingBox(), actor1,
            freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceComputeNormals3D<T, Descriptor>, lattice2.getBoundingBox(), actor2,
            freeSurfaceArgs2, pl);

        /***** New level ******/
        pl++;

        if (useSurfaceTension1) {
            integrateProcessingFunctional(
                new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngle1),
                lattice1.getBoundingBox(), actor1, freeSurfaceArgs1, pl);
        }
        if (useSurfaceTension2) {
            integrateProcessingFunctional(
                new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngle2),
                lattice2.getBoundingBox(), actor2, freeSurfaceArgs2, pl);
        }

        integrateProcessingFunctional(
            new FreeSurfaceMassChange3D<T, Descriptor>, lattice1.getBoundingBox(), actor1,
            freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceMassChange3D<T, Descriptor>, lattice2.getBoundingBox(), actor2,
            freeSurfaceArgs2, pl);

        integrateProcessingFunctional(
            new FreeSurfaceCompletion3D<T, Descriptor>, lattice1.getBoundingBox(), actor1,
            freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceCompletion3D<T, Descriptor>, lattice2.getBoundingBox(), actor2,
            freeSurfaceArgs2, pl);

        integrateProcessingFunctional(
            new FreeSurfaceMacroscopic3D<T, Descriptor>(incompressibleModel1),
            lattice1.getBoundingBox(), actor1, freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceMacroscopic3D<T, Descriptor>(incompressibleModel2),
            lattice2.getBoundingBox(), actor2, freeSurfaceArgs2, pl);

        /***** New level ******/
        pl++;

        bool useRepellingForceCoupling = !util::isZero(interactionStrength);

        if (useRepellingForceCoupling) {
            // integrateProcessingFunctional (
            //     new MultiFreeSurfaceRepellingForceCoupling3D<T,Descriptor>(interactionStrength,
            //     rhoDefault1, rhoDefault2), lattice2.getBoundingBox(), actor2,
            //     multiFreeSurfaceArgs, pl );
            integrateProcessingFunctional(
                new MultiFreeSurfaceComplexCoupling3D<T, Descriptor>(
                    interactionStrength, rhoDefault1, rhoDefault2),
                lattice2.getBoundingBox(), actor2, multiFreeSurfaceArgs, pl);
        } else {
            integrateProcessingFunctional(
                new MultiFreeSurfaceVelocityContinuityCoupling3D<T, Descriptor>(
                    rhoDefault1, rhoDefault2),
                lattice2.getBoundingBox(), actor2, multiFreeSurfaceArgs, pl);
        }

        // integrateProcessingFunctional (
        //     new MultiFreeSurfaceOneWayCoupling3D<T,Descriptor>(interactionStrength,rhoDefault1,
        //     rhoDefault2), lattice2.getBoundingBox(), actor2, multiFreeSurfaceArgs, pl );

        /***** New level ******/
        if (useSurfaceTension1 || useSurfaceTension2)
            pl++;

        if (useSurfaceTension1) {
            integrateProcessingFunctional(
                new FreeSurfaceAddSurfaceTension3D<T, Descriptor>(
                    surfaceTension1, incompressibleModel1),
                lattice1.getBoundingBox(), actor1, freeSurfaceArgs1, pl);
        }
        if (useSurfaceTension2) {
            integrateProcessingFunctional(
                new FreeSurfaceAddSurfaceTension3D<T, Descriptor>(
                    surfaceTension2, incompressibleModel2),
                lattice2.getBoundingBox(), actor2, freeSurfaceArgs2, pl);
        }

        integrateProcessingFunctional(
            new FreeSurfaceStabilize3D<T, Descriptor>(), lattice1.getBoundingBox(), actor1,
            freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceStabilize3D<T, Descriptor>(), lattice2.getBoundingBox(), actor2,
            freeSurfaceArgs2, pl);

        /***** New level ******/  // Maybe no new level is necessary here ...
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), lattice1.getBoundingBox(),
            actor1, freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), lattice2.getBoundingBox(),
            actor2, freeSurfaceArgs2, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoDefault1),
            lattice1.getBoundingBox(), actor1, freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoDefault2),
            lattice2.getBoundingBox(), actor2, freeSurfaceArgs2, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(dynamics1->clone(), force1),
            lattice1.getBoundingBox(), actor1, freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(dynamics2->clone(), force2),
            lattice2.getBoundingBox(), actor2, freeSurfaceArgs2, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoDefault1),
            lattice1.getBoundingBox(), actor1, freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoDefault2),
            lattice2.getBoundingBox(), actor2, freeSurfaceArgs2, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>(),
            lattice1.getBoundingBox(), actor1, freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>(),
            lattice2.getBoundingBox(), actor2, freeSurfaceArgs2, pl);

        integrateProcessingFunctional(
            new FreeSurfaceComputeStatistics3D<T, Descriptor>, lattice1.getBoundingBox(), actor1,
            freeSurfaceArgs1, pl);
        integrateProcessingFunctional(
            new FreeSurfaceComputeStatistics3D<T, Descriptor>, lattice2.getBoundingBox(), actor2,
            freeSurfaceArgs2, pl);

        bool useForce1 = !util::isZero(norm(force1));
        if (useForce1) {
            integrateProcessingFunctional(
                new FreeSurfaceAddExternalForce3D<T, Descriptor>(rhoDefault1),
                lattice1.getBoundingBox(), actor1, freeSurfaceArgs1, pl);
        }
        bool useForce2 = !util::isZero(norm(force2));
        if (useForce2) {
            integrateProcessingFunctional(
                new FreeSurfaceAddExternalForce3D<T, Descriptor>(rhoDefault2),
                lattice2.getBoundingBox(), actor2, freeSurfaceArgs2, pl);
        }
    }

    Dynamics<T, Descriptor> *dynamics1, *dynamics2;
    T rhoDefault1, rhoDefault2;
    T surfaceTension1, surfaceTension2;
    T contactAngle1, contactAngle2;
    bool useSurfaceTension1, useSurfaceTension2;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force1, force2;
    MultiBlockLattice3D<T, Descriptor> lattice1, lattice2;
    MultiContainerBlock3D helperLists1, helperLists2;
    MultiScalarField3D<T> mass1, mass2;
    MultiScalarField3D<int> flag1, flag2;
    MultiScalarField3D<T> volumeFraction1, volumeFraction2;
    MultiScalarField3D<T> curvature1, curvature2;
    MultiScalarField3D<T> outsideDensity1, outsideDensity2;
    MultiScalarField3D<T> rhoBar1, rhoBar2;
    MultiTensorField3D<T, 3> j1, j2;
    MultiTensorField3D<T, 3> normal1, normal2;
    T interactionStrength;
    std::vector<MultiBlock3D *> rhoBarJparam1, rhoBarJparam2;
    std::vector<MultiBlock3D *> freeSurfaceArgs1, freeSurfaceArgs2, multiFreeSurfaceArgs;
};

}  // namespace plb

#endif  // MULTI_FREE_SURFACE_MODEL_3D_H
