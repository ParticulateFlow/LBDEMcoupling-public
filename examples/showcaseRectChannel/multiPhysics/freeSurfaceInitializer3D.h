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

#ifndef FREE_SURFACE_INITIALIZER_3D_H
#define FREE_SURFACE_INITIALIZER_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiPhysics/freeSurfaceUtil3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
class DefaultInitializeFreeSurface3D : public BoxProcessingFunctional3D {
public:
    DefaultInitializeFreeSurface3D(
        Dynamics<T, Descriptor> *dynamicsTemplate_,
        Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_, T rhoIni_ = (T)1.,
        bool useRhoIni_ = true, bool useZeroMomentum_ = true, bool initializeCell_ = true) :
        dynamicsTemplate(dynamicsTemplate_),
        emptyNodeDynamicsTemplate(emptyNodeDynamicsTemplate_),
        force(force_),
        rhoIni(rhoIni_),
        useRhoIni(useRhoIni_),
        useZeroMomentum(useZeroMomentum_),
        initializeCell(initializeCell_)
    { }
    DefaultInitializeFreeSurface3D(
        Dynamics<T, Descriptor> *dynamicsTemplate_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_, T rhoIni_ = (T)1.,
        bool useRhoIni_ = true, bool useZeroMomentum_ = true, bool initializeCell_ = true) :
        dynamicsTemplate(dynamicsTemplate_),
        force(force_),
        rhoIni(rhoIni_),
        useRhoIni(useRhoIni_),
        useZeroMomentum(useZeroMomentum_),
        initializeCell(initializeCell_)
    {
        emptyNodeDynamicsTemplate = new NoDynamics<T, Descriptor>(rhoIni);
    }
    DefaultInitializeFreeSurface3D(DefaultInitializeFreeSurface3D<T, Descriptor> const &rhs) :
        dynamicsTemplate(rhs.dynamicsTemplate->clone()),
        emptyNodeDynamicsTemplate(rhs.emptyNodeDynamicsTemplate->clone()),
        force(rhs.force),
        rhoIni(rhs.rhoIni),
        useRhoIni(rhs.useRhoIni),
        useZeroMomentum(rhs.useZeroMomentum),
        initializeCell(rhs.initializeCell)
    { }
    DefaultInitializeFreeSurface3D<T, Descriptor> *operator=(
        DefaultInitializeFreeSurface3D<T, Descriptor> const &rhs)
    {
        DefaultInitializeFreeSurface3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(DefaultInitializeFreeSurface3D<T, Descriptor> &rhs)
    {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(emptyNodeDynamicsTemplate, rhs.emptyNodeDynamicsTemplate);
        std::swap(force, rhs.force);
        std::swap(rhoIni, rhs.rhoIni);
        std::swap(useRhoIni, rhs.useRhoIni);
        std::swap(useZeroMomentum, rhs.useZeroMomentum);
        std::swap(initializeCell, rhs.initializeCell);
    }
    virtual ~DefaultInitializeFreeSurface3D()
    {
        delete dynamicsTemplate;
        delete emptyNodeDynamicsTemplate;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual DefaultInitializeFreeSurface3D<T, Descriptor> *clone() const
    {
        return new DefaultInitializeFreeSurface3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::dataStructure;    // Fluid
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass
        modified[4] = modif::staticVariables;  // Volume-fraction
        modified[5] = modif::staticVariables;  // Flag-status
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::staticVariables;  // Outside density.
    }

private:
    Dynamics<T, Descriptor> *dynamicsTemplate;
    Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force;
    T rhoIni;
    bool useRhoIni, useZeroMomentum, initializeCell;
};

// Same as DefaultInitializeFreeSurface, but without initializing the Volume-fraction.
template <typename T, template <typename U> class Descriptor>
class PartiallyDefaultInitializeFreeSurface3D : public BoxProcessingFunctional3D {
public:
    PartiallyDefaultInitializeFreeSurface3D(
        Dynamics<T, Descriptor> *dynamicsTemplate_,
        Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_, T rhoIni_ = (T)1.,
        bool useRhoIni_ = true, bool useZeroMomentum_ = true, bool initializeCell_ = true) :
        dynamicsTemplate(dynamicsTemplate_),
        emptyNodeDynamicsTemplate(emptyNodeDynamicsTemplate_),
        force(force_),
        rhoIni(rhoIni_),
        useRhoIni(useRhoIni_),
        useZeroMomentum(useZeroMomentum_),
        initializeCell(initializeCell_)
    { }
    PartiallyDefaultInitializeFreeSurface3D(
        Dynamics<T, Descriptor> *dynamicsTemplate_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_, T rhoIni_ = (T)1.,
        bool useRhoIni_ = true, bool useZeroMomentum_ = true, bool initializeCell_ = true) :
        dynamicsTemplate(dynamicsTemplate_),
        force(force_),
        rhoIni(rhoIni_),
        useRhoIni(useRhoIni_),
        useZeroMomentum(useZeroMomentum_),
        initializeCell(initializeCell_)
    {
        emptyNodeDynamicsTemplate = new NoDynamics<T, Descriptor>(rhoIni);
    }
    PartiallyDefaultInitializeFreeSurface3D(
        PartiallyDefaultInitializeFreeSurface3D<T, Descriptor> const &rhs) :
        dynamicsTemplate(rhs.dynamicsTemplate->clone()),
        emptyNodeDynamicsTemplate(rhs.emptyNodeDynamicsTemplate->clone()),
        force(rhs.force),
        rhoIni(rhs.rhoIni),
        useRhoIni(rhs.useRhoIni),
        useZeroMomentum(rhs.useZeroMomentum),
        initializeCell(rhs.initializeCell)
    { }
    PartiallyDefaultInitializeFreeSurface3D<T, Descriptor> *operator=(
        PartiallyDefaultInitializeFreeSurface3D<T, Descriptor> const &rhs)
    {
        PartiallyDefaultInitializeFreeSurface3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(PartiallyDefaultInitializeFreeSurface3D<T, Descriptor> &rhs)
    {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(emptyNodeDynamicsTemplate, rhs.emptyNodeDynamicsTemplate);
        std::swap(force, rhs.force);
        std::swap(rhoIni, rhs.rhoIni);
        std::swap(useRhoIni, rhs.useRhoIni);
        std::swap(useZeroMomentum, rhs.useZeroMomentum);
        std::swap(initializeCell, rhs.initializeCell);
    }
    virtual ~PartiallyDefaultInitializeFreeSurface3D()
    {
        delete dynamicsTemplate;
        delete emptyNodeDynamicsTemplate;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual PartiallyDefaultInitializeFreeSurface3D<T, Descriptor> *clone() const
    {
        return new PartiallyDefaultInitializeFreeSurface3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::dataStructure;    // Fluid
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass
        modified[4] = modif::staticVariables;  // Volume-fraction
        modified[5] = modif::staticVariables;  // Flag-status
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::staticVariables;  // Outside density.
    }

private:
    Dynamics<T, Descriptor> *dynamicsTemplate;
    Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force;
    T rhoIni;
    bool useRhoIni, useZeroMomentum, initializeCell;
};

template <typename T, class InsideFunction>
class AnalyticalIniVolumeFraction3D : public BoxProcessingFunctional3D {
public:
    AnalyticalIniVolumeFraction3D(InsideFunction const &insideFunction_, plint subDivision_ = 5) :
        insideFunction(insideFunction_), subDivision(subDivision_)
    {
        PLB_ASSERT(subDivision > 1);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual AnalyticalIniVolumeFraction3D<T, InsideFunction> *clone() const
    {
        return new AnalyticalIniVolumeFraction3D<T, InsideFunction>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;  // Volume-fraction
        modified[1] = modif::staticVariables;  // Flag-status
    }

private:
    void subDomainVolumeFraction(plint iX, plint iY, plint iZ, int &flag, T &volumeFraction);

private:
    InsideFunction const &insideFunction;
    plint subDivision;
};

template <typename T, class InsideFunction>
void analyticalIniVolumeFraction(
    MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flagStatus,
    InsideFunction const &insideFunction, Box3D domain, plint subDivision = 5);

template <typename T, class InsideFunction>
void analyticalIniVolumeFraction(
    MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flagStatus,
    InsideFunction const &insideFunction, plint subDivision = 5);

template <typename T, template <typename U> class Descriptor>
class ConstantIniVelocityFreeSurface3D : public BoxProcessingFunctional3D {
public:
    ConstantIniVelocityFreeSurface3D(Array<T, 3> velocity_, T rhoIni_) :
        velocity(velocity_), rhoIni(rhoIni_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual ConstantIniVelocityFreeSurface3D<T, Descriptor> *clone() const
    {
        return new ConstantIniVelocityFreeSurface3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables;  // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::nothing;          // Volume-fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    Array<T, 3> velocity;
    T rhoIni;
};

template <typename T, template <typename U> class Descriptor>
class InletConstVolumeFraction3D : public BoxProcessingFunctional3D {
public:
    InletConstVolumeFraction3D(T volumeFraction_) : volumeFraction(volumeFraction_) { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual InletConstVolumeFraction3D<T, Descriptor> *clone() const
    {
        return new InletConstVolumeFraction3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume-fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    T volumeFraction;
};

template <typename T, template <typename U> class Descriptor>
class MaskedInletConstVolumeFraction3D : public BoxProcessingFunctional3D {
public:
    MaskedInletConstVolumeFraction3D(T volumeFraction_, int whichFlag_) :
        volumeFraction(volumeFraction_), whichFlag(whichFlag_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual MaskedInletConstVolumeFraction3D<T, Descriptor> *clone() const
    {
        return new MaskedInletConstVolumeFraction3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;          // rhoBar.
        modified[1] = modif::staticVariables;  // mass.
        modified[2] = modif::nothing;          // mask.
    }

private:
    T volumeFraction;
    int whichFlag;
};

template <typename T, template <typename U> class Descriptor>
void maskedInletConstVolumeFraction3D(
    MultiScalarField3D<T> &rhoBar, MultiScalarField3D<T> &mass, MultiScalarField3D<int> &mask,
    T volumeFraction, int whichFlag, Box3D domain);

template <typename T, template <typename U> class Descriptor>
class N_MaskedInletConstVolumeFraction3D : public BoxProcessingFunctional3D {
public:
    N_MaskedInletConstVolumeFraction3D(T volumeFraction_) : volumeFraction(volumeFraction_) { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual N_MaskedInletConstVolumeFraction3D<T, Descriptor> *clone() const
    {
        return new N_MaskedInletConstVolumeFraction3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;          // rhoBar.
        modified[1] = modif::staticVariables;  // mass.
        modified[2] = modif::nothing;          // mask.
    }

private:
    T volumeFraction;
};

template <typename T, template <typename U> class Descriptor>
class OutletMaximumVolumeFraction3D : public BoxProcessingFunctional3D {
public:
    OutletMaximumVolumeFraction3D(T volumeFraction_) : volumeFraction(volumeFraction_) { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual OutletMaximumVolumeFraction3D<T, Descriptor> *clone() const
    {
        return new OutletMaximumVolumeFraction3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume-fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    T volumeFraction;
};

template <typename T, template <typename U> class Descriptor>
class N_MaskedOutletMaximumVolumeFraction3D : public BoxProcessingFunctional3D {
public:
    N_MaskedOutletMaximumVolumeFraction3D(T volumeFraction_) : volumeFraction(volumeFraction_) { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual N_MaskedOutletMaximumVolumeFraction3D<T, Descriptor> *clone() const
    {
        return new N_MaskedOutletMaximumVolumeFraction3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;          // rhoBar.
        modified[1] = modif::staticVariables;  // mass.
        modified[2] = modif::nothing;          // mask.
    }

private:
    T volumeFraction;
};

template <typename T, template <typename U> class Descriptor>
class OutletVolumeFractionInRange3D : public BoxProcessingFunctional3D {
public:
    OutletVolumeFractionInRange3D(T minFraction_, T maxFraction_) :
        minFraction(minFraction_), maxFraction(maxFraction_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual OutletVolumeFractionInRange3D<T, Descriptor> *clone() const
    {
        return new OutletVolumeFractionInRange3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume-fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    T minFraction, maxFraction;
};

template <typename T, template <typename U> class Descriptor>
class OutletMaximumVolumeFraction2_3D : public BoxProcessingFunctional3D {
public:
    OutletMaximumVolumeFraction2_3D(T volumeFraction_) : volumeFraction(volumeFraction_) { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual OutletMaximumVolumeFraction2_3D<T, Descriptor> *clone() const
    {
        return new OutletMaximumVolumeFraction2_3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables;  // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume-fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    T volumeFraction;
};

template <typename T, template <typename U> class Descriptor>
class NoSlipMaximumVolumeFraction3D : public BoxProcessingFunctional3D {
public:
    NoSlipMaximumVolumeFraction3D(T volumeFraction_) : volumeFraction(volumeFraction_) { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual NoSlipMaximumVolumeFraction3D<T, Descriptor> *clone() const
    {
        return new NoSlipMaximumVolumeFraction3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::nothing;          // Volume-fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    T volumeFraction;
};

template <typename T, template <typename U> class Descriptor>
class PunchSphere3D : public BoxProcessingFunctional3D {
public:
    PunchSphere3D(
        Array<T, 3> const &center_, T radius_, T rho0_,
        Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate_) :
        center(center_),
        radius(radius_),
        rho0(rho0_),
        emptyNodeDynamicsTemplate(emptyNodeDynamicsTemplate_)
    { }
    PunchSphere3D(Array<T, 3> const &center_, T radius_, T rho0_) :
        center(center_), radius(radius_), rho0(rho0_)
    {
        emptyNodeDynamicsTemplate = new NoDynamics<T, Descriptor>(rho0);
    }
    PunchSphere3D(PunchSphere3D<T, Descriptor> const &rhs) :
        center(rhs.center),
        radius(rhs.radius),
        rho0(rhs.rho0),
        emptyNodeDynamicsTemplate(rhs.emptyNodeDynamicsTemplate->clone())
    { }
    PunchSphere3D<T, Descriptor> *operator=(PunchSphere3D<T, Descriptor> const &rhs)
    {
        PunchSphere3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(PunchSphere3D<T, Descriptor> &rhs)
    {
        std::swap(center, rhs.center);
        std::swap(radius, rhs.radius);
        std::swap(rho0, rhs.rho0);
        std::swap(emptyNodeDynamicsTemplate, rhs.emptyNodeDynamicsTemplate);
    }
    virtual ~PunchSphere3D()
    {
        delete emptyNodeDynamicsTemplate;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual PunchSphere3D<T, Descriptor> *clone() const
    {
        return new PunchSphere3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::dataStructure;    // Fluid.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::staticVariables;  // Volume-fraction.
        modified[5] = modif::staticVariables;  // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::staticVariables;  // Outside density.
    }

private:
    Array<T, 3> center;
    T radius;
    T rho0;
    Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate;
};

template <typename T, template <typename U> class Descriptor>
class AnalyticalPunchSphere3D : public BoxProcessingFunctional3D {
public:
    AnalyticalPunchSphere3D(
        Array<T, 3> const &center_, T radius_, T rho0_,
        Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate_, plint subDivision_ = 5) :
        center(center_),
        radius(radius_),
        rho0(rho0_),
        subDivision(subDivision_),
        emptyNodeDynamicsTemplate(emptyNodeDynamicsTemplate_)
    {
        PLB_ASSERT(subDivision > 1);
    }
    AnalyticalPunchSphere3D(
        Array<T, 3> const &center_, T radius_, T rho0_, plint subDivision_ = 5) :
        center(center_), radius(radius_), rho0(rho0_), subDivision(subDivision_)
    {
        PLB_ASSERT(subDivision > 1);
        emptyNodeDynamicsTemplate = new NoDynamics<T, Descriptor>(rho0);
    }
    AnalyticalPunchSphere3D(AnalyticalPunchSphere3D<T, Descriptor> const &rhs) :
        center(rhs.center),
        radius(rhs.radius),
        rho0(rhs.rho0),
        subDivision(rhs.subDivision),
        emptyNodeDynamicsTemplate(rhs.emptyNodeDynamicsTemplate->clone())
    { }
    AnalyticalPunchSphere3D<T, Descriptor> *operator=(
        AnalyticalPunchSphere3D<T, Descriptor> const &rhs)
    {
        AnalyticalPunchSphere3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(AnalyticalPunchSphere3D<T, Descriptor> &rhs)
    {
        std::swap(center, rhs.center);
        std::swap(radius, rhs.radius);
        std::swap(rho0, rhs.rho0);
        std::swap(subDivision, rhs.subDivision);
        std::swap(emptyNodeDynamicsTemplate, rhs.emptyNodeDynamicsTemplate);
    }
    virtual ~AnalyticalPunchSphere3D()
    {
        delete emptyNodeDynamicsTemplate;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual AnalyticalPunchSphere3D<T, Descriptor> *clone() const
    {
        return new AnalyticalPunchSphere3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::dataStructure;    // Fluid.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::staticVariables;  // Volume-fraction.
        modified[5] = modif::staticVariables;  // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::staticVariables;  // Outside density.
    }

private:
    bool isInsideSphere(T x, T y, T z)
    {
        return (normSqr(Array<T, 3>(x, y, z) - center) < radius * radius);
    }
    void subDomainVolumeFraction(
        plint globalX, plint globalY, plint globalZ, int &flag, T &volumeFraction);

private:
    Array<T, 3> center;
    T radius;
    T rho0;
    plint subDivision;
    Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate;
};

template <typename T, template <typename U> class Descriptor>
class CalculateAverageSphereDensity3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    CalculateAverageSphereDensity3D(Array<T, 3> const &center_, T radius_) :
        center(center_), radius(radius_), averageDensityId(this->getStatistics().subscribeAverage())
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual CalculateAverageSphereDensity3D<T, Descriptor> *clone() const
    {
        return new CalculateAverageSphereDensity3D<T, Descriptor>(*this);
    }
    T getAverageDensity() const
    {
        return this->getStatistics().getAverage(averageDensityId);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
    }

private:
    Array<T, 3> center;
    T radius;
    plint averageDensityId;
};

}  // namespace plb

#endif  // FREE_SURFACE_INITIALIZER_3D_H
