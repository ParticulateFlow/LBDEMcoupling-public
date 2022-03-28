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

#ifndef TWO_PHASE_MODEL_3D_H
#define TWO_PHASE_MODEL_3D_H

#include <memory>

#include "core/globalDefs.h"
#include "core/util.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiPhysics/freeSurfaceModel3D.h"

namespace plb {

template <typename T, int nDim>
class IniFilteredDensity3D : public BoxProcessingFunctional3D_ST<T, T, nDim> {
public:
    virtual void process(
        Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField);
    virtual IniFilteredDensity3D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, int nDim>
void IniFilteredDensity3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (int iDim = 0; iDim < nDim; ++iDim) {
                    tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iDim] =
                        scalarField.get(iX, iY, iZ);
                }
            }
        }
    }
}

template <typename T, int nDim>
IniFilteredDensity3D<T, nDim> *IniFilteredDensity3D<T, nDim>::clone() const
{
    return new IniFilteredDensity3D<T, nDim>(*this);
}

template <typename T, int nDim>
void IniFilteredDensity3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
void iniFilteredDensity(
    MultiTensorField3D<T, nDim> &densities, MultiScalarField3D<T> &density, Box3D domain)
{
    applyProcessingFunctional(new IniFilteredDensity3D<T, nDim>(), domain, density, densities);
}

template <typename T, int nDim>
class AddFilteredDensity3D : public BoxProcessingFunctional3D_ST<T, T, nDim> {
public:
    virtual void process(
        Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField);
    virtual AddFilteredDensity3D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, int nDim>
void AddFilteredDensity3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                for (int iDim = 0; iDim < nDim - 1; ++iDim) {
                    tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iDim] =
                        tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iDim + 1];
                }
                tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z)[nDim - 1] =
                    scalarField.get(iX, iY, iZ);
            }
        }
    }
}

template <typename T, int nDim>
AddFilteredDensity3D<T, nDim> *AddFilteredDensity3D<T, nDim>::clone() const
{
    return new AddFilteredDensity3D<T, nDim>(*this);
}

template <typename T, int nDim>
void AddFilteredDensity3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, int nDim>
void addFilteredDensity(
    MultiTensorField3D<T, nDim> &densities, MultiScalarField3D<T> &density, Box3D domain)
{
    applyProcessingFunctional(new AddFilteredDensity3D<T, nDim>(), domain, density, densities);
}

template <typename T, int nDim>
class GetFilteredDensity3D : public BoxProcessingFunctional3D_ST<T, T, nDim> {
public:
    virtual void process(
        Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField);
    virtual GetFilteredDensity3D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, int nDim>
void GetFilteredDensity3D<T, nDim>::process(
    Box3D domain, ScalarField3D<T> &scalarField, TensorField3D<T, nDim> &tensorField)
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T averageDensity = T();
                for (int iDim = 0; iDim < nDim; ++iDim) {
                    averageDensity +=
                        tensorField.get(iX + offset.x, iY + offset.y, iZ + offset.z)[iDim];
                }
                scalarField.get(iX, iY, iZ) = averageDensity;
            }
        }
    }
}

template <typename T, int nDim>
GetFilteredDensity3D<T, nDim> *GetFilteredDensity3D<T, nDim>::clone() const
{
    return new GetFilteredDensity3D<T, nDim>(*this);
}

template <typename T, int nDim>
void GetFilteredDensity3D<T, nDim>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template <typename T, int nDim>
void getFilteredDensity(
    MultiTensorField3D<T, nDim> &densities, MultiScalarField3D<T> &density, Box3D domain)
{
    applyProcessingFunctional(new GetFilteredDensity3D<T, nDim>(), domain, density, densities);
}

/** kinetic: The pressure on the interface is entirely
 *               determined by the mass exchange between either phase.
 *  dynamic:     The pressure is considered to be a fluctuation around
 *               "outsideDensity", which can be either equal to rhoEmpty, or
 *               to the result of the pressure-correction model.
 *  constRho:    The pressure is constant throughout the bubble interface
 *               and takes its value from "outsideDensity", which can be either
 *               equal to rhoEmpty, or to the result of the pressure-correction model.
 **/
typedef enum {
    undefined = -1,
    kinetic = 1,
    dynamic = 2,
    bubblePressure = 3,
    constRho = 4,
    freeSurface = 5
} TwoPhaseModel;

inline TwoPhaseModel stringToTwoPhaseModel(std::string modelName)
{
    if (modelName == "kinetic") {
        return kinetic;
    } else if (modelName == "dynamic") {
        return dynamic;
    } else if (modelName == "bubblePressure") {
        return bubblePressure;
    } else if (modelName == "constRho") {
        return constRho;
    } else if (modelName == "freeSurface") {
        return freeSurface;
    } else {
        PLB_ASSERT(false);
    }

    return undefined;
}

template <typename T, template <typename U> class Descriptor>
class TwoPhaseAveragePressure3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    TwoPhaseAveragePressure3D(T densityRatio_, T rhoDefault_, TwoPhaseModel model_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseAveragePressure3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseAveragePressure3D<T, Descriptor>(*this);
    }
    T getAveragePressure() const;
    T getAveragePressure2() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::nothing;          // Fluid 2.
            modified[11] = modif::nothing;          // rhoBar2.
            modified[12] = modif::nothing;          // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    plint sumRho1_ID, sumRho2_ID, weight1_ID, weight2_ID;
    T densityRatio, rhoDefault;
    TwoPhaseModel model;
};

template <typename T, template <typename U> class Descriptor>
class TwoPhaseAverageVelocity3D : public PlainReductiveBoxProcessingFunctional3D {
public:
    TwoPhaseAverageVelocity3D(TwoPhaseModel model_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseAverageVelocity3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseAverageVelocity3D<T, Descriptor>(*this);
    }
    Array<T, 3> getAverageVelocity() const;
    Array<T, 3> getAverageVelocity2() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::nothing;          // Fluid 2.
            modified[11] = modif::nothing;          // rhoBar2.
            modified[12] = modif::nothing;          // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    Array<plint, 3> sumVel1_ID, sumVel2_ID;
    plint weight1_ID, weight2_ID;
    TwoPhaseModel model;
};

template <typename T, template <typename U> class Descriptor>
class TwoPhaseComputePressure3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseComputePressure3D(T densityRatio_, T rhoDefault_, TwoPhaseModel model_) :
        densityRatio(densityRatio_),
        rhoDefault(rhoDefault_),
        model(model_),
        computeFluid1(true),
        computeFluid2(true)
    { }
    TwoPhaseComputePressure3D(
        T densityRatio_, T rhoDefault_, TwoPhaseModel model_, bool computeFluid1_,
        bool computeFluid2_) :
        densityRatio(densityRatio_),
        rhoDefault(rhoDefault_),
        model(model_),
        computeFluid1(computeFluid1_),
        computeFluid2(computeFluid2_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseComputePressure3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseComputePressure3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;  // Fluid.
        modified[1] = modif::nothing;  // rhoBar.
        modified[2] = modif::nothing;  // j.
        modified[3] = modif::nothing;  // Mass.
        modified[4] = modif::nothing;  // Volume fraction.
        modified[5] = modif::nothing;  // Flag-status.
        modified[6] = modif::nothing;  // Normal.
        modified[7] = modif::nothing;  // Interface lists.
        modified[8] = modif::nothing;  // Curvature.
        modified[9] = modif::nothing;  // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::nothing;          // Fluid 2.
            modified[11] = modif::nothing;          // rhoBar2.
            modified[12] = modif::nothing;          // j2.
            modified[13] = modif::nothing;          // mass2.
            modified[14] = modif::staticVariables;  // resulting density
        } else {
            modified[10] = modif::staticVariables;  // resulting density
        }
    }

private:
    T densityRatio, rhoDefault;
    TwoPhaseModel model;
    bool computeFluid1;
    bool computeFluid2;
};

template <typename T, template <typename U> class Descriptor>
class TwoPhaseComputeVelocity3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseComputeVelocity3D(T densityRatio_, bool useFreeSurfaceLimit_) :
        densityRatio(densityRatio_),
        computeFluid1(true),
        computeFluid2(true),
        useFreeSurfaceLimit(useFreeSurfaceLimit_)
    { }
    TwoPhaseComputeVelocity3D(
        T densityRatio_, bool computeFluid1_, bool computeFluid2_, bool useFreeSurfaceLimit_) :
        densityRatio(densityRatio_),
        computeFluid1(computeFluid1_),
        computeFluid2(computeFluid2_),
        useFreeSurfaceLimit(useFreeSurfaceLimit_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseComputeVelocity3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseComputeVelocity3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;  // Fluid.
        modified[1] = modif::nothing;  // rhoBar.
        modified[2] = modif::nothing;  // j.
        modified[3] = modif::nothing;  // Mass.
        modified[4] = modif::nothing;  // Volume fraction.
        modified[5] = modif::nothing;  // Flag-status.
        modified[6] = modif::nothing;  // Normal.
        modified[7] = modif::nothing;  // Interface lists.
        modified[8] = modif::nothing;  // Curvature.
        modified[9] = modif::nothing;  // Outside density.
        if (!useFreeSurfaceLimit) {
            modified[10] = modif::nothing;          // Fluid 2.
            modified[11] = modif::nothing;          // rhoBar2.
            modified[12] = modif::nothing;          // j2.
            modified[13] = modif::nothing;          // mass2.
            modified[14] = modif::staticVariables;  // resulting velocity
        } else {
            modified[10] = modif::staticVariables;  // resulting velocity
        }
    }

private:
    T densityRatio;
    bool computeFluid1;
    bool computeFluid2;
    bool useFreeSurfaceLimit;
};

template <typename T, template <typename U> class Descriptor>
class TwoPhasePunchSphere3D : public BoxProcessingFunctional3D {
public:
    TwoPhasePunchSphere3D(
        Array<T, 3> const &center_, T radius_, T rhoEmpty_, T referenceDensity_, T densityRatio_,
        Dynamics<T, Descriptor> *dynamicsTemplate2_, TwoPhaseModel model_) :
        center(center_),
        radius(radius_),
        rhoEmpty(rhoEmpty_),
        referenceDensity(referenceDensity_),
        densityRatio(densityRatio_),
        dynamicsTemplate2(dynamicsTemplate2_),
        model(model_)
    { }
    TwoPhasePunchSphere3D(TwoPhasePunchSphere3D<T, Descriptor> const &rhs) :
        center(rhs.center),
        radius(rhs.radius),
        rhoEmpty(rhs.rhoEmpty),
        referenceDensity(rhs.referenceDensity),
        densityRatio(rhs.densityRatio),
        dynamicsTemplate2(rhs.dynamicsTemplate2->clone()),
        model(rhs.model)
    { }
    TwoPhasePunchSphere3D<T, Descriptor> &operator=(TwoPhasePunchSphere3D<T, Descriptor> const &rhs)
    {
        TwoPhasePunchSphere3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(TwoPhasePunchSphere3D<T, Descriptor> &rhs)
    {
        std::swap(center, rhs.center);
        std::swap(radius, rhs.radius);
        std::swap(rhoEmpty, rhs.rhoEmpty);
        std::swap(referenceDensity, rhs.referenceDensity);
        std::swap(densityRatio, rhs.densityRatio);
        std::swap(dynamicsTemplate2, rhs.dynamicsTemplate2);
        std::swap(model, rhs.model);
    }
    ~TwoPhasePunchSphere3D()
    {
        delete dynamicsTemplate2;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhasePunchSphere3D<T, Descriptor> *clone() const
    {
        return new TwoPhasePunchSphere3D<T, Descriptor>(*this);
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
        if (model != freeSurface) {
            modified[10] = modif::dataStructure;    // Fluid 2.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    Array<T, 3> center;
    T radius;
    T rhoEmpty, referenceDensity, densityRatio;
    Dynamics<T, Descriptor> *dynamicsTemplate2;
    TwoPhaseModel model;
};

template <typename T, template <typename U> class Descriptor>
class TwoPhasePunchRectangle3D : public BoxProcessingFunctional3D {
public:
    TwoPhasePunchRectangle3D(
        Box3D rectangle_, T rhoEmpty_, T referenceDensity_, T densityRatio_,
        Dynamics<T, Descriptor> *dynamicsTemplate2_, TwoPhaseModel model_) :
        rectangle(rectangle_),
        rhoEmpty(rhoEmpty_),
        referenceDensity(referenceDensity_),
        densityRatio(densityRatio_),
        dynamicsTemplate2(dynamicsTemplate2_),
        model(model_)
    { }
    TwoPhasePunchRectangle3D(TwoPhasePunchRectangle3D<T, Descriptor> const &rhs) :
        rectangle(rhs.rectangle),
        rhoEmpty(rhs.rhoEmpty),
        referenceDensity(rhs.referenceDensity),
        densityRatio(rhs.densityRatio),
        dynamicsTemplate2(rhs.dynamicsTemplate2->clone()),
        model(rhs.model)
    { }
    TwoPhasePunchRectangle3D<T, Descriptor> &operator=(
        TwoPhasePunchRectangle3D<T, Descriptor> const &rhs)
    {
        TwoPhasePunchRectangle3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(TwoPhasePunchRectangle3D<T, Descriptor> &rhs)
    {
        std::swap(rectangle, rhs.rectangle);
        std::swap(rhoEmpty, rhs.rhoEmpty);
        std::swap(referenceDensity, rhs.referenceDensity);
        std::swap(densityRatio, rhs.densityRatio);
        std::swap(dynamicsTemplate2, rhs.dynamicsTemplate2);
        std::swap(model, rhs.model);
    }
    ~TwoPhasePunchRectangle3D()
    {
        delete dynamicsTemplate2;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhasePunchRectangle3D<T, Descriptor> *clone() const
    {
        return new TwoPhasePunchRectangle3D<T, Descriptor>(*this);
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
        if (model != freeSurface) {
            modified[10] = modif::dataStructure;    // Fluid 2.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    Box3D rectangle;
    T rhoEmpty, referenceDensity, densityRatio;
    Dynamics<T, Descriptor> *dynamicsTemplate2;
    TwoPhaseModel model;
};

template <typename T, template <typename U> class Descriptor>
void twoPhasePunchSphere(
    std::vector<MultiBlock3D *> const &twoPhaseArgs, Array<T, 3> const &center, T radius,
    T rhoEmpty, T referenceDensity, T densityRatio, Dynamics<T, Descriptor> &dynamics,
    TwoPhaseModel model, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void twoPhasePunchRectangle(
    std::vector<MultiBlock3D *> const &twoPhaseArgs, Box3D rectangle, T rhoEmpty,
    T referenceDensity, T densityRatio, Dynamics<T, Descriptor> &dynamics, TwoPhaseModel model,
    Box3D domain);

template <typename T, template <typename U> class Descriptor>
T computeAverageSphereDensity(
    std::vector<MultiBlock3D *> const &twoPhaseArgs, Array<T, 3> const &center, T radius,
    Box3D domain);

/// Data structure for holding lists of cells along the free surface in an AtomicContainerBlock.
template <typename T, template <typename U> class Descriptor>
struct TwoPhaseInterfaceLists : public ContainerBlockData {
    typedef Array<plint, Descriptor<T>::d> Node;
    struct ExtrapolInfo {
        ExtrapolInfo() : density(T())
        {
            j.resetToZero();
            PiNeq.resetToZero();
        }
        T density;
        Array<T, 3> j;
        Array<T, SymmetricTensor<T, Descriptor>::n> PiNeq;
    };
    /// Holds all nodes which have excess mass.
    std::map<Node, T> massExcess;
    /// Holds all nodes which have excess mass for fluid 2.
    std::map<Node, T> massExcess2;
    /// Holds all nodes that need to change status from interface to fluid.
    std::set<Node> interfaceToFluid;
    /// Holds all nodes that need to change status from interface to empty.
    std::set<Node> interfaceToEmpty;
    /// Holds all nodes that need to change status from empty to interface.
    std::map<Node, ExtrapolInfo> emptyToInterface;
    /// Holds all nodes that need to change status from fluid to interface.
    std::map<Node, ExtrapolInfo> fluidToInterface;

    virtual TwoPhaseInterfaceLists<T, Descriptor> *clone() const
    {
        return new TwoPhaseInterfaceLists<T, Descriptor>(*this);
    }
};

/// A wrapper offering convenient access to the free-surface data provided to
/// data processors. Avoids verbous casting, asserting, etc.
template <typename T, template <typename U> class Descriptor>
class TwoPhaseProcessorParam3D {
public:
    typedef typename TwoPhaseInterfaceLists<T, Descriptor>::Node Node;
    typedef typename TwoPhaseInterfaceLists<T, Descriptor>::ExtrapolInfo ExtrapolInfo;
    TwoPhaseProcessorParam3D(std::vector<AtomicBlock3D *> &atomicBlocks)
    {
        if (atomicBlocks.size() >= 14) {
            useFreeSurfaceLimit = false;
        } else {
            useFreeSurfaceLimit = true;
        }

        fluid_ = dynamic_cast<BlockLattice3D<T, Descriptor> *>(atomicBlocks[0]);
        PLB_ASSERT(fluid_);

        rhoBar_ = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[1]);
        PLB_ASSERT(rhoBar_);

        j_ = dynamic_cast<TensorField3D<T, 3> *>(atomicBlocks[2]);
        PLB_ASSERT(j_);

        mass_ = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[3]);
        PLB_ASSERT(mass_);

        volumeFraction_ = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[4]);
        PLB_ASSERT(volumeFraction_);

        flag_ = dynamic_cast<ScalarField3D<int> *>(atomicBlocks[5]);
        PLB_ASSERT(flag_);

        normal_ = dynamic_cast<TensorField3D<T, 3> *>(atomicBlocks[6]);
        PLB_ASSERT(normal_);

        containerInterfaceLists_ = dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[7]);
        PLB_ASSERT(containerInterfaceLists_);

        interfaceLists_ = dynamic_cast<TwoPhaseInterfaceLists<T, Descriptor> *>(
            containerInterfaceLists_->getData());
        // PLB_ASSERT(interfaceLists_);
        // Put the assertion at the usage of interfaceLists, so we can still work with both
        // freeSurfaceProcessorParam and twoPhaseProcessorParam.

        curvature_ = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[8]);
        PLB_ASSERT(curvature_);

        outsideDensity_ = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[9]);
        PLB_ASSERT(outsideDensity_);

        if (!useFreeSurfaceLimit) {
            fluid2_ = dynamic_cast<BlockLattice3D<T, Descriptor> *>(atomicBlocks[10]);
            PLB_ASSERT(fluid2_);

            rhoBar2_ = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[11]);
            PLB_ASSERT(rhoBar2_);

            j2_ = dynamic_cast<TensorField3D<T, 3> *>(atomicBlocks[12]);
            PLB_ASSERT(j2_);

            mass2_ = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[13]);
            PLB_ASSERT(mass2_);
        } else {
            fluid2_ = 0;
            rhoBar2_ = 0;
            j2_ = 0;
            mass2_ = 0;
        }

        absoluteOffset = fluid_->getLocation();
        relativeOffsetRhoBar = computeRelativeDisplacement(*fluid_, *rhoBar_);
        relativeOffsetJ = computeRelativeDisplacement(*fluid_, *j_);
        relativeOffsetMass = computeRelativeDisplacement(*fluid_, *mass_);
        relativeOffsetVF = computeRelativeDisplacement(*fluid_, *volumeFraction_);
        relativeOffsetFS = computeRelativeDisplacement(*fluid_, *flag_);
        relativeOffsetNormal = computeRelativeDisplacement(*fluid_, *normal_);
        relativeOffsetC = computeRelativeDisplacement(*fluid_, *curvature_);
        relativeOffsetOD = computeRelativeDisplacement(*fluid_, *outsideDensity_);

        if (!useFreeSurfaceLimit) {
            relativeOffsetFluid2 = computeRelativeDisplacement(*fluid_, *fluid2_);
            relativeOffsetRhoBar2 = computeRelativeDisplacement(*fluid_, *rhoBar2_);
            relativeOffsetJ2 = computeRelativeDisplacement(*fluid_, *j2_);
            relativeOffsetMass2 = computeRelativeDisplacement(*fluid_, *mass2_);
        }
    }
    Cell<T, Descriptor> &cell(plint iX, plint iY, plint iZ)
    {
        return fluid_->get(iX, iY, iZ);
    }
    Cell<T, Descriptor> &cell2(plint iX, plint iY, plint iZ)
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        return fluid2_->get(
            iX + relativeOffsetFluid2.x, iY + relativeOffsetFluid2.y, iZ + relativeOffsetFluid2.z);
    }
    T &mass(plint iX, plint iY, plint iZ)
    {
        return mass_->get(
            iX + relativeOffsetMass.x, iY + relativeOffsetMass.y, iZ + relativeOffsetMass.z);
    }
    T &mass2(plint iX, plint iY, plint iZ)
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        return mass2_->get(
            iX + relativeOffsetMass2.x, iY + relativeOffsetMass2.y, iZ + relativeOffsetMass2.z);
    }
    T &volumeFraction(plint iX, plint iY, plint iZ)
    {
        return volumeFraction_->get(
            iX + relativeOffsetVF.x, iY + relativeOffsetVF.y, iZ + relativeOffsetVF.z);
    }
    T &curvature(plint iX, plint iY, plint iZ)
    {
        return curvature_->get(
            iX + relativeOffsetC.x, iY + relativeOffsetC.y, iZ + relativeOffsetC.z);
    }
    T &outsideDensity(plint iX, plint iY, plint iZ)
    {
        return outsideDensity_->get(
            iX + relativeOffsetOD.x, iY + relativeOffsetOD.y, iZ + relativeOffsetOD.z);
    }
    int &flag(plint iX, plint iY, plint iZ)
    {
        return flag_->get(
            iX + relativeOffsetFS.x, iY + relativeOffsetFS.y, iZ + relativeOffsetFS.z);
    }
    void setForce(
        plint iX, plint iY, plint iZ,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> const &force)
    {
        setExternalForce(cell(iX, iY, iZ), force);
    }
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> getForce(plint iX, plint iY, plint iZ)
    {
        return getExternalForce(cell(iX, iY, iZ));
    }
    void setForce2(
        plint iX, plint iY, plint iZ,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> const &force)
    {
        setExternalForce(cell2(iX, iY, iZ), force);
    }
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> getForce2(plint iX, plint iY, plint iZ)
    {
        return getExternalForce(cell2(iX, iY, iZ));
    }
    void setMomentum(plint iX, plint iY, plint iZ, Array<T, 3> const &momentum)
    {
        j_->get(iX + relativeOffsetJ.x, iY + relativeOffsetJ.y, iZ + relativeOffsetJ.z) = momentum;
    }
    void setMomentum2(plint iX, plint iY, plint iZ, Array<T, 3> const &momentum)
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        j2_->get(iX + relativeOffsetJ2.x, iY + relativeOffsetJ2.y, iZ + relativeOffsetJ2.z) =
            momentum;
    }
    Array<T, 3> getMomentum(plint iX, plint iY, plint iZ)
    {
        return j_->get(iX + relativeOffsetJ.x, iY + relativeOffsetJ.y, iZ + relativeOffsetJ.z);
    }
    Array<T, 3> getMomentum2(plint iX, plint iY, plint iZ)
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        return j2_->get(iX + relativeOffsetJ2.x, iY + relativeOffsetJ2.y, iZ + relativeOffsetJ2.z);
    }
    T getDensity(plint iX, plint iY, plint iZ)
    {
        return Descriptor<T>::fullRho(rhoBar_->get(
            iX + relativeOffsetRhoBar.x, iY + relativeOffsetRhoBar.y, iZ + relativeOffsetRhoBar.z));
    }
    T getDensity2(plint iX, plint iY, plint iZ)
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        return Descriptor<T>::fullRho(rhoBar2_->get(
            iX + relativeOffsetRhoBar2.x, iY + relativeOffsetRhoBar2.y,
            iZ + relativeOffsetRhoBar2.z));
    }
    void setDensity(plint iX, plint iY, plint iZ, T rho)
    {
        rhoBar_->get(
            iX + relativeOffsetRhoBar.x, iY + relativeOffsetRhoBar.y, iZ + relativeOffsetRhoBar.z) =
            Descriptor<T>::rhoBar(rho);
    }
    void setDensity2(plint iX, plint iY, plint iZ, T rho)
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        rhoBar2_->get(
            iX + relativeOffsetRhoBar2.x, iY + relativeOffsetRhoBar2.y,
            iZ + relativeOffsetRhoBar2.z) = Descriptor<T>::rhoBar(rho);
    }
    void setNormal(plint iX, plint iY, plint iZ, Array<T, 3> const &normal)
    {
        normal_->get(
            iX + relativeOffsetNormal.x, iY + relativeOffsetNormal.y, iZ + relativeOffsetNormal.z) =
            normal;
    }
    Array<T, 3> getNormal(plint iX, plint iY, plint iZ)
    {
        return normal_->get(
            iX + relativeOffsetNormal.x, iY + relativeOffsetNormal.y, iZ + relativeOffsetNormal.z);
    }

    void attributeDynamics(plint iX, plint iY, plint iZ, Dynamics<T, Descriptor> *dynamics)
    {
        fluid_->attributeDynamics(iX, iY, iZ, dynamics);
    }

    void attributeDynamics2(plint iX, plint iY, plint iZ, Dynamics<T, Descriptor> *dynamics)
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        fluid2_->attributeDynamics(iX, iY, iZ, dynamics);
    }

    bool isBoundary(plint iX, plint iY, plint iZ)
    {
        return cell(iX, iY, iZ).getDynamics().isBoundary();
    }

    void addToTotalMass(T addedTotalMass)
    {
        fluid_->getInternalStatistics().gatherSum(0, addedTotalMass);
    }
    void addToLostMass(T addedLostMass)
    {
        fluid_->getInternalStatistics().gatherSum(1, addedLostMass);
    }
    void addToTotalMass2(T addedTotalMass2)
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        fluid2_->getInternalStatistics().gatherSum(0, addedTotalMass2);
    }
    void addToLostMass2(T addedLostMass2)
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        fluid2_->getInternalStatistics().gatherSum(1, addedLostMass2);
    }
    void addToInterfaceCells(plint addedInterfaceCells)
    {
        fluid_->getInternalStatistics().gatherIntSum(0, addedInterfaceCells);
    }
    T getSumMassMatrix() const
    {
        return fluid_->getInternalStatistics().getSum(0);
    }
    T getSumLostMass() const
    {
        return fluid_->getInternalStatistics().getSum(1);
    }
    T getTotalMass() const
    {
        return getSumMassMatrix() + getSumLostMass();
    }
    T getSumMassMatrix2() const
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        return fluid2_->getInternalStatistics().getSum(0);
    }
    T getSumLostMass2() const
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        return fluid2_->getInternalStatistics().getSum(1);
    }
    T getTotalMass2() const
    {
        PLB_ASSERT(!useFreeSurfaceLimit);
        return getSumMassMatrix2() + getSumLostMass2();
    }
    plint getNumInterfaceCells() const
    {
        return fluid_->getInternalStatistics().getIntSum(0);
    }

    std::map<Node, T> &massExcess()
    {
        PLB_ASSERT(interfaceLists_);
        return interfaceLists_->massExcess;
    }
    std::map<Node, T> &massExcess2()
    {
        PLB_ASSERT(interfaceLists_);
        return interfaceLists_->massExcess2;
    }
    std::set<Node> &interfaceToFluid()
    {
        PLB_ASSERT(interfaceLists_);
        return interfaceLists_->interfaceToFluid;
    }
    std::set<Node> &interfaceToEmpty()
    {
        PLB_ASSERT(interfaceLists_);
        return interfaceLists_->interfaceToEmpty;
    }
    std::map<Node, ExtrapolInfo> &emptyToInterface()
    {
        PLB_ASSERT(interfaceLists_);
        return interfaceLists_->emptyToInterface;
    }
    std::map<Node, ExtrapolInfo> &fluidToInterface()
    {
        PLB_ASSERT(interfaceLists_);
        return interfaceLists_->fluidToInterface;
    }

    Dot3D const &absOffset() const
    {
        return absoluteOffset;
    }
    Box3D getBoundingBox() const
    {
        return volumeFraction_->getBoundingBox();
    }

private:
    bool useFreeSurfaceLimit;
    BlockLattice3D<T, Descriptor> *fluid_, *fluid2_;
    ScalarField3D<T> *rhoBar_, *rhoBar2_;
    TensorField3D<T, 3> *j_, *j2_;
    ScalarField3D<T> *mass_, *mass2_;
    ScalarField3D<T> *volumeFraction_;
    ScalarField3D<int> *flag_;
    TensorField3D<T, 3> *normal_;
    AtomicContainerBlock3D *containerInterfaceLists_;
    TwoPhaseInterfaceLists<T, Descriptor> *interfaceLists_;
    ScalarField3D<T> *curvature_;
    ScalarField3D<T> *outsideDensity_;

    Dot3D absoluteOffset, relativeOffsetFluid2, relativeOffsetRhoBar, relativeOffsetRhoBar2,
        relativeOffsetJ, relativeOffsetJ2, relativeOffsetMass, relativeOffsetMass2,
        relativeOffsetVF, relativeOffsetFS, relativeOffsetNormal, relativeOffsetC, relativeOffsetOD;
};

/// Create a parameter-list for most free-surface data processors.
template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlock3D *> aggregateTwoPhaseParams(
    MultiBlockLattice3D<T, Descriptor> &fluid, MultiBlockLattice3D<T, Descriptor> *fluid2,
    MultiScalarField3D<T> &rhoBar, MultiScalarField3D<T> *rhoBar2, MultiTensorField3D<T, 3> &j,
    MultiTensorField3D<T, 3> *j2, MultiScalarField3D<T> &mass, MultiScalarField3D<T> *mass2,
    MultiScalarField3D<T> &volumeFraction, MultiScalarField3D<int> &flag,
    MultiTensorField3D<T, 3> &normal, MultiContainerBlock3D &interfaceLists,
    MultiScalarField3D<T> &curvature, MultiScalarField3D<T> &outsideDensity)
{
    std::vector<MultiBlock3D *> aggregation;

    aggregation.push_back(&fluid);
    aggregation.push_back(&rhoBar);
    aggregation.push_back(&j);
    aggregation.push_back(&mass);
    aggregation.push_back(&volumeFraction);
    aggregation.push_back(&flag);
    aggregation.push_back(&normal);
    aggregation.push_back(&interfaceLists);
    aggregation.push_back(&curvature);
    aggregation.push_back(&outsideDensity);
    if (fluid2) {
        PLB_ASSERT(rhoBar2);
        PLB_ASSERT(j2);
        PLB_ASSERT(mass2);
        aggregation.push_back(fluid2);
        aggregation.push_back(rhoBar2);
        aggregation.push_back(j2);
        aggregation.push_back(mass2);
    }

    return aggregation;
}

template <typename T, template <typename U> class Descriptor>
class DefaultInitializeTwoPhase3D : public BoxProcessingFunctional3D {
public:
    DefaultInitializeTwoPhase3D(
        Dynamics<T, Descriptor> *dynamicsTemplate_, Dynamics<T, Descriptor> *dynamicsTemplate2_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force2_, T rhoIni_,
        bool useFreeSurfaceLimit_) :
        dynamicsTemplate(dynamicsTemplate_),
        dynamicsTemplate2(dynamicsTemplate2_),
        force(force_),
        force2(force2_),
        rhoIni(rhoIni_),
        useFreeSurfaceLimit(useFreeSurfaceLimit_)
    { }
    DefaultInitializeTwoPhase3D(DefaultInitializeTwoPhase3D<T, Descriptor> const &rhs) :
        dynamicsTemplate(rhs.dynamicsTemplate->clone()),
        dynamicsTemplate2(rhs.dynamicsTemplate2->clone()),
        force(rhs.force),
        force2(rhs.force2),
        rhoIni(rhs.rhoIni),
        useFreeSurfaceLimit(rhs.useFreeSurfaceLimit)
    { }
    DefaultInitializeTwoPhase3D<T, Descriptor> &operator=(
        DefaultInitializeTwoPhase3D<T, Descriptor> const &rhs)
    {
        DefaultInitializeTwoPhase3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(DefaultInitializeTwoPhase3D<T, Descriptor> &rhs)
    {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(dynamicsTemplate2, rhs.dynamicsTemplate2);
        std::swap(force, rhs.force);
        std::swap(force2, rhs.force2);
        std::swap(rhoIni, rhs.rhoIni);
        std::swap(useFreeSurfaceLimit, rhs.useFreeSurfaceLimit);
    }
    virtual ~DefaultInitializeTwoPhase3D()
    {
        delete dynamicsTemplate;
        delete dynamicsTemplate2;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual DefaultInitializeTwoPhase3D<T, Descriptor> *clone() const
    {
        return new DefaultInitializeTwoPhase3D(*this);
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
        modified[9] = modif::nothing;          // Outside density.
        if (!useFreeSurfaceLimit) {
            modified[10] = modif::dataStructure;    // Fluid 2.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    Dynamics<T, Descriptor> *dynamicsTemplate, *dynamicsTemplate2;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force, force2;
    T rhoIni;
    bool useFreeSurfaceLimit;
};

template <typename T, template <typename U> class Descriptor>
class PartiallyDefaultInitializeTwoPhase3D : public BoxProcessingFunctional3D {
public:
    PartiallyDefaultInitializeTwoPhase3D(
        Dynamics<T, Descriptor> *dynamicsTemplate_, Dynamics<T, Descriptor> *dynamicsTemplate2_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force2_, T rhoIni_,
        bool useFreeSurfaceLimit_) :
        dynamicsTemplate(dynamicsTemplate_),
        dynamicsTemplate2(dynamicsTemplate2_),
        force(force_),
        force2(force2_),
        rhoIni(rhoIni_),
        useFreeSurfaceLimit(useFreeSurfaceLimit_)
    { }
    PartiallyDefaultInitializeTwoPhase3D(
        PartiallyDefaultInitializeTwoPhase3D<T, Descriptor> const &rhs) :
        dynamicsTemplate(rhs.dynamicsTemplate->clone()),
        dynamicsTemplate2(rhs.dynamicsTemplate2->clone()),
        force(rhs.force),
        force2(rhs.force2),
        rhoIni(rhs.rhoIni),
        useFreeSurfaceLimit(rhs.useFreeSurfaceLimit)
    { }
    PartiallyDefaultInitializeTwoPhase3D<T, Descriptor> &operator=(
        PartiallyDefaultInitializeTwoPhase3D<T, Descriptor> const &rhs)
    {
        PartiallyDefaultInitializeTwoPhase3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(PartiallyDefaultInitializeTwoPhase3D<T, Descriptor> &rhs)
    {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(dynamicsTemplate2, rhs.dynamicsTemplate2);
        std::swap(force, rhs.force);
        std::swap(force2, rhs.force2);
        std::swap(rhoIni, rhs.rhoIni);
        std::swap(useFreeSurfaceLimit, rhs.useFreeSurfaceLimit);
    }
    virtual ~PartiallyDefaultInitializeTwoPhase3D()
    {
        delete dynamicsTemplate;
        delete dynamicsTemplate2;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual PartiallyDefaultInitializeTwoPhase3D<T, Descriptor> *clone() const
    {
        return new PartiallyDefaultInitializeTwoPhase3D(*this);
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
        modified[9] = modif::nothing;          // Outside density.
        if (!useFreeSurfaceLimit) {
            modified[10] = modif::dataStructure;    // Fluid 2.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    Dynamics<T, Descriptor> *dynamicsTemplate, *dynamicsTemplate2;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force, force2;
    T rhoIni;
    bool useFreeSurfaceLimit;
};

// Functional to impose an initial constant velocity to one or to both fluids.
// CAUTION: This data processor must be called only after proper initialization of the
// TwoPhaseFields3D data structure.
template <typename T, template <typename U> class Descriptor, class Function>
class ConstantIniVelocityTwoPhase3D : public BoxProcessingFunctional3D {
public:
    ConstantIniVelocityTwoPhase3D(
        Array<T, 3> velocity_, Function f_, bool imposeToFluid1_, bool imposeToFluid2_,
        bool useFreeSurfaceLimit_) :
        velocity(velocity_),
        f(f_),
        imposeToFluid1(imposeToFluid1_),
        imposeToFluid2(imposeToFluid2_),
        useFreeSurfaceLimit(useFreeSurfaceLimit_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual ConstantIniVelocityTwoPhase3D<T, Descriptor, Function> *clone() const
    {
        return new ConstantIniVelocityTwoPhase3D<T, Descriptor, Function>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables;  // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if (!useFreeSurfaceLimit) {
            modified[10] = modif::staticVariables;  // Fluid 2.
            modified[11] = modif::nothing;          // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::nothing;          // mass2.
        }
    }

private:
    Array<T, 3> velocity;
    Function f;
    bool imposeToFluid1, imposeToFluid2;
    bool useFreeSurfaceLimit;
};

/// Compute the mass balance on every node in the domain, and store in mass matrix.
/** Input:
 *   - Flag-status:   needed in bulk+1
 *   - Mass:          needed in bulk
 *   - Volume fraction: needed in bulk
 *   - Populations:   needed in bulk+1
 * Output:
 *   - mass.
 **/
template <typename T, template <typename U> class Descriptor>
class TwoPhaseMassChange3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseMassChange3D(bool useFreeSurfaceLimit_) : useFreeSurfaceLimit(useFreeSurfaceLimit_) { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseMassChange3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseMassChange3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if (!useFreeSurfaceLimit) {
            modified[10] = modif::nothing;          // Fluid 2.
            modified[11] = modif::nothing;          // rhoBar2.
            modified[12] = modif::nothing;          // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    bool useFreeSurfaceLimit;
};

/// Completion scheme on the post-collide populations on interface cells.
/** Input:
 *   - Flag-status:   needed in bulk+1
 *   - Volume fraction: needed in bulk+1
 *   - Populations:   needed in bulk+1
 *   - Momentum:      needed in bulk+1
 *   - Density:       needed in bulk+1
 * Output:
 *   - Populations.
 **/
// ASK: This data processor loops over the whole volume. Is this really
//      necessary, or could one of the lists be used instead?
template <typename T, template <typename U> class Descriptor>
class TwoPhaseCompletion3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseCompletion3D(bool useFreeSurfaceLimit_) : useFreeSurfaceLimit(useFreeSurfaceLimit_) { }
    virtual TwoPhaseCompletion3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseCompletion3D<T, Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;  // Fluid. Should be: staticVariables.
        modified[1] = modif::nothing;  // rhoBar.
        modified[2] = modif::nothing;  // j.
        modified[3] = modif::nothing;  // Mass.
        modified[4] = modif::nothing;  // Volume fraction.
        modified[5] = modif::nothing;  // Flag-status.
        modified[6] = modif::nothing;  // Normal.
        modified[7] = modif::nothing;  // Interface-lists.
        modified[8] = modif::nothing;  // Curvature.
        modified[9] = modif::nothing;  // Outside density.
        if (!useFreeSurfaceLimit) {
            modified[10] = modif::nothing;  // Fluid 2. Should be: staticVariables.
            modified[11] = modif::nothing;  // rhoBar2.
            modified[12] = modif::nothing;  // j2.
            modified[13] = modif::nothing;  // mass2.
        }
    }

private:
    bool useFreeSurfaceLimit;
};

/// Compute and store mass-fraction and macroscopic variables.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Mass:          needed in bulk
 *   - Populations:   needed in bulk
 * Output:
 *   - mass-fraction, density, momentum, flag (because setting bounce-back).
 **/
template <typename T, template <typename U> class Descriptor>
class TwoPhaseMacroscopic3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseMacroscopic3D(T rhoDefault_, T densityRatio_, T surfaceTension_, TwoPhaseModel model_) :
        rhoDefault(rhoDefault_),
        densityRatio(densityRatio_),
        surfaceTension(surfaceTension_),
        model(model_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseMacroscopic3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseMacroscopic3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables;  // Fluid. Should be: staticVariables.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass. Should be: staticVariables.
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::staticVariables;  // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::staticVariables;  // Fluid 2. Should be: staticVariables.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    T rhoDefault;
    T densityRatio;
    T surfaceTension;
    TwoPhaseModel model;
};

/// Stabilization scheme on the post-collide populations on interface cells.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Volume fraction: needed in bulk
 *   - Populations:   needed in bulk
 *   - Momentum:      needed in bulk
 *   - Density:       needed in bulk
 * Output:
 *   - Populations.
 *   - Momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class TwoPhaseStabilize3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseStabilize3D(TwoPhaseModel model_) : model(model_) { }
    virtual TwoPhaseStabilize3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseStabilize3D<T, Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables;  // Fluid. Should be: staticVariables.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::staticVariables;  // Fluid 2. Should be: staticVariables.
            modified[11] = modif::nothing;          // rhoBar2.
            modified[12] = modif::nothing;          // j2.
            modified[13] = modif::nothing;          // mass2.
        }
    }

private:
    T rhoDefault;
    TwoPhaseModel model;
};

template <typename T, template <typename U> class Descriptor>
class VerifyTwoPhase : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual VerifyTwoPhase<T, Descriptor> *clone() const
    {
        return new VerifyTwoPhase<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;   // Fluid. Should be: nothing.
        modified[1] = modif::nothing;   // rhoBar.
        modified[2] = modif::nothing;   // j.
        modified[3] = modif::nothing;   // Mass. Should be: nothing.
        modified[4] = modif::nothing;   // Volume fraction.
        modified[5] = modif::nothing;   // Flag-status.
        modified[6] = modif::nothing;   // Normal.
        modified[7] = modif::nothing;   // Interface-lists.
        modified[8] = modif::nothing;   // Curvature.
        modified[9] = modif::nothing;   // Outside density.
        modified[10] = modif::nothing;  // Fluid 2. Should be: nothing.
        modified[11] = modif::nothing;  // rhoBar2.
        modified[12] = modif::nothing;  // j2.
        modified[13] = modif::nothing;  // mass2.
    }
};

template <typename T, template <typename U> class Descriptor>
class TwoPhaseInterfaceFilter : public BoxProcessingFunctional3D {
public:
    TwoPhaseInterfaceFilter(T rhoDefault_, T densityRatio_, TwoPhaseModel model_) :
        rhoDefault(rhoDefault_), densityRatio(densityRatio_), model(model_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseInterfaceFilter<T, Descriptor> *clone() const
    {
        return new TwoPhaseInterfaceFilter<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables;   // Fluid. Should be: staticVariables.
        modified[1] = modif::staticVariables;   // rhoBar.
        modified[2] = modif::staticVariables;   // j.
        modified[3] = modif::staticVariables;   // Mass. Should be: staticVariables.
        modified[4] = modif::staticVariables;   // Volume fraction.
        modified[5] = modif::staticVariables;   // Flag-status.
        modified[6] = modif::staticVariables;   // Normal.
        modified[7] = modif::staticVariables;   // Interface-lists.
        modified[8] = modif::staticVariables;   // Curvature.
        modified[9] = modif::staticVariables;   // Outside density.
        modified[10] = modif::dataStructure;    // Fluid 2. Should be: staticVariables.
        modified[11] = modif::staticVariables;  // rhoBar2.
        modified[12] = modif::staticVariables;  // j2.
        modified[13] = modif::staticVariables;  // mass2.
    }

private:
    T rhoDefault;
    T densityRatio;
    TwoPhaseModel model;
};

/** Input:
 *   - interface-to-fluid list: needed in bulk+1
 *   - interface-to-empty list: needed in bulk+1
 *   - density: needed in bulk+1
 *   - mass:    needed in bulk+1
 *   - flag:    needed in bulk+1
 * Output:
 *   - flag, dynamics, mass, volumeFraction, density, force, momentum
 *   - mass-excess-list: defined in bulk+1
 **/
template <typename T, template <typename U> class Descriptor>
class TwoPhaseIniInterfaceToAnyNodes3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseIniInterfaceToAnyNodes3D(
        T rhoDefault_, Dynamics<T, Descriptor> *dynamicsTemplate_,
        Dynamics<T, Descriptor> *dynamicsTemplate2_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> const &force_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> const &force2_, T densityRatio_,
        T surfaceTension_, TwoPhaseModel model_) :
        rhoDefault(rhoDefault_),
        dynamicsTemplate(dynamicsTemplate_),
        dynamicsTemplate2(dynamicsTemplate2_),
        force(force_),
        force2(force2_),
        densityRatio(densityRatio_),
        surfaceTension(surfaceTension_),
        model(model_)
    { }
    TwoPhaseIniInterfaceToAnyNodes3D(TwoPhaseIniInterfaceToAnyNodes3D<T, Descriptor> const &rhs) :
        rhoDefault(rhs.rhoDefault),
        dynamicsTemplate(rhs.dynamicsTemplate->clone()),
        dynamicsTemplate2(rhs.dynamicsTemplate2->clone()),
        force(rhs.force),
        force2(rhs.force2),
        densityRatio(rhs.densityRatio),
        surfaceTension(rhs.surfaceTension),
        model(rhs.model)
    { }
    TwoPhaseIniInterfaceToAnyNodes3D<T, Descriptor> &operator=(
        TwoPhaseIniInterfaceToAnyNodes3D<T, Descriptor> const &rhs)
    {
        TwoPhaseIniInterfaceToAnyNodes3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(TwoPhaseIniInterfaceToAnyNodes3D<T, Descriptor> &rhs)
    {
        std::swap(rhoDefault, rhs.rhoDefault);
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(dynamicsTemplate2, rhs.dynamicsTemplate2);
        std::swap(force, rhs.force);
        std::swap(force2, rhs.force2);
        std::swap(densityRatio, rhs.densityRatio);
        std::swap(surfaceTension, rhs.surfaceTension);
        std::swap(model, rhs.model);
    }
    virtual ~TwoPhaseIniInterfaceToAnyNodes3D()
    {
        delete dynamicsTemplate;
        delete dynamicsTemplate2;
    }
    TwoPhaseIniInterfaceToAnyNodes3D(T rhoDefault_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);

    virtual TwoPhaseIniInterfaceToAnyNodes3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseIniInterfaceToAnyNodes3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] =
            modif::nothing;  // Fluid. Gets assigned new dynamics. Should be: dataStructure
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::nothing;          // j. Should be: staticVariables.
        modified[3] = modif::staticVariables;  // Mass. Is redistributed and initialized from
                                               // neighborying density.
        modified[4] =
            modif::nothing;  // Volume fraction. Is default-initialized. Should be: staticVariables.
        modified[5] =
            modif::staticVariables;    // Flag-status. Is adapted according to cell-change lists.
        modified[6] = modif::nothing;  // Normal.
        modified[7] = modif::nothing;  // Interface-lists. Read-only.
        modified[8] = modif::nothing;  // Curvature.
        modified[9] = modif::nothing;  // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::nothing;          // Fluid 2. Should be: dataStructure
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::nothing;          // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    T rhoDefault;
    Dynamics<T, Descriptor> *dynamicsTemplate, *dynamicsTemplate2;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force, force2;
    T densityRatio, surfaceTension;
    TwoPhaseModel model;
};

/// Based on the previously computed empty->interface list, initialize flow variables for
///   new interface cells.
/** Input:
 *   - Populations: needed in bulk+0
 *   - Momentum:    needed in bulk+1
 *   - Density:     needed in bulk+1
 *   - Flag-status: needed in bulk+0
 * Output:
 *   - flag-status:   initialized to "interface" on corresponding cells.
 *   - lattice:       initialized from neighbor averages on new interface cells.
 *   - mass:          initialized to zero on new interface cells.
 *   - mass-fraction: initialized to zero on new interface cells.
 *   - momentum
 **/
template <typename T, template <typename U> class Descriptor>
class TwoPhaseIniEmptyToInterfaceNodes3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseIniEmptyToInterfaceNodes3D(
        Dynamics<T, Descriptor> *dynamicsTemplate_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_, T densityRatio_,
        TwoPhaseModel model_) :
        dynamicsTemplate(dynamicsTemplate_),
        force(force_),
        densityRatio(densityRatio_),
        model(model_)
    { }
    TwoPhaseIniEmptyToInterfaceNodes3D(
        TwoPhaseIniEmptyToInterfaceNodes3D<T, Descriptor> const &rhs) :
        dynamicsTemplate(rhs.dynamicsTemplate->clone()),
        force(rhs.force),
        densityRatio(rhs.densityRatio),
        model(rhs.model)
    { }
    TwoPhaseIniEmptyToInterfaceNodes3D<T, Descriptor> &operator=(
        TwoPhaseIniEmptyToInterfaceNodes3D<T, Descriptor> const &rhs)
    {
        TwoPhaseIniEmptyToInterfaceNodes3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(TwoPhaseIniEmptyToInterfaceNodes3D<T, Descriptor> &rhs)
    {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(force, rhs.force);
        std::swap(densityRatio, rhs.densityRatio);
        std::swap(model, rhs.model);
    }
    ~TwoPhaseIniEmptyToInterfaceNodes3D()
    {
        delete dynamicsTemplate;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseIniEmptyToInterfaceNodes3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseIniEmptyToInterfaceNodes3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid. Should be: dataStructure
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::nothing;          // j. Should be: staticVariables.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;  // Volume fraction, read-only. Should be: staticVariables
        modified[5] = modif::staticVariables;  // Flag-status, read-only.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;  // Interface-lists. Read access to gasCellToInitializeData.
        modified[8] = modif::nothing;  // Curvature.
        modified[9] = modif::nothing;  // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::dataStructure;    // Fluid 2.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::nothing;          // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    Dynamics<T, Descriptor> *dynamicsTemplate;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce>
        force;  // Body force, for initialization of the new interface cell.
    T densityRatio;
    TwoPhaseModel model;
};

/// Isolated cells cannot be part of the interface. This data processor spots and
/// removes them.
/** Input:
 *   - Flag-status: needed in bulk+2
 *   - mass:        needed in bulk+1
 *   - density:     needed in bulk+1
 * Output:
 *   - interfaceToFluidNodes:   initialized in bulk+1
 *   - interfaceToEmptyNodes:   initialized in bulk+1
 *   - massExcess list:         initialized in bulk+1
 *   - mass, density, mass-fraction, dynamics, force, momentum, flag: in bulk+1
 **/
template <typename T, template <typename U> class Descriptor>
class TwoPhaseRemoveFalseInterfaceCells3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseRemoveFalseInterfaceCells3D(T rhoDefault_, TwoPhaseModel model_) :
        rhoDefault(rhoDefault_), model(model_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseRemoveFalseInterfaceCells3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseRemoveFalseInterfaceCells3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;  // Fluid: Gets NoDynamics when node changes to empty. Should
                                       // be: dataStructure.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::nothing;          // j. Should be: staticVariables.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume fraction. Should be: staticVariables.
        modified[5] = modif::staticVariables;  // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::nothing;          // Fluid 2.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::nothing;          // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    T rhoDefault;
    TwoPhaseModel model;
};

template <typename T, template <typename U> class Descriptor>
class TwoPhaseInitializeInterfaceLists3D : public BoxProcessingFunctional3D {
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
    {
        PLB_ASSERT(atomicBlocks.size() == 1);

        AtomicContainerBlock3D *containerInterfaceLists =
            dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[0]);
        PLB_ASSERT(containerInterfaceLists);
        TwoPhaseInterfaceLists<T, Descriptor> *interfaceLists =
            new TwoPhaseInterfaceLists<T, Descriptor>;
        containerInterfaceLists->setData(interfaceLists);
    }
    virtual TwoPhaseInitializeInterfaceLists3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseInitializeInterfaceLists3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        // Default-assign potential other parameters present in a multi-fluid system.
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables;
    }
};

/// Wrapper for execution of TwoPhaseInitializeInterfaceLists3D.
template <typename T, template <typename U> class Descriptor>
void twoPhaseInitializeInterfaceLists3D(MultiContainerBlock3D &interfaceListBlock)
{
    std::vector<MultiBlock3D *> arg;
    arg.push_back(&interfaceListBlock);
    applyProcessingFunctional(
        new TwoPhaseInitializeInterfaceLists3D<T, Descriptor>, interfaceListBlock.getBoundingBox(),
        arg);
}

template <typename T, template <typename U> class Descriptor>
class TwoPhaseComputeInterfaceLists3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseComputeInterfaceLists3D(
        TwoPhaseModel model_, T rhoDefault_, T densityRatio_, T surfaceTension_) :
        model(model_),
        rhoDefault(rhoDefault_),
        densityRatio(densityRatio_),
        surfaceTension(surfaceTension_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseComputeInterfaceLists3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseComputeInterfaceLists3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;  // Fluid.
        modified[1] = modif::nothing;  // rhoBar.
        modified[2] = modif::nothing;  // j.
        modified[3] = modif::nothing;  // Mass.
        modified[4] = modif::nothing;  // Volume fraction.
        modified[5] = modif::nothing;  // Flag-status.
        modified[6] = modif::nothing;  // Normal.
        modified[7] = modif::nothing;  // Interface-lists.
        modified[8] = modif::nothing;  // Curvature.
        modified[9] = modif::nothing;  // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::nothing;  // Fluid 2.
            modified[11] = modif::nothing;  // rhoBar2.
            modified[12] = modif::nothing;  // j2.
            modified[13] = modif::nothing;  // mass2.
        }
    }

private:
    static T kappa;  // Safety threshold for state-change, to prevent back-and-forth oscillations.
    TwoPhaseModel model;
    T rhoDefault;
    T densityRatio;
    T surfaceTension;
};

/// Enforce exact mass balance when interface cells become fluid or empty.
/** Input:
 *   - mass-excess list: needed in bulk+1
 *   - Flag-status: needed in bulk+2
 *   - mass:        needed in bulk+2
 *   - density:     needed in bulk+2
 * Output:
 *   - mass, mass-fraction
 **/
template <typename T, template <typename U> class Descriptor>
class TwoPhaseEqualMassExcessReDistribution3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseEqualMassExcessReDistribution3D(TwoPhaseModel model_) : model(model_) { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseEqualMassExcessReDistribution3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseEqualMassExcessReDistribution3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::dataStructure;    // Fluid.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::dataStructure;    // Fluid 2.
            modified[11] = modif::staticVariables;  // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::staticVariables;  // mass2.
        }
    }

private:
    TwoPhaseModel model;
};

/// Addition of external forces.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Volume fraction: needed in bulk
 *   - Populations:   needed in bulk
 *   - Momentum:      needed in bulk
 *   - Density:       needed in bulk
 * Output:
 *   - Momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class TwoPhaseAddExternalForce3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseAddExternalForce3D(T rhoDefault_, TwoPhaseModel model_) :
        rhoDefault(rhoDefault_), model(model_)
    { }
    virtual TwoPhaseAddExternalForce3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseAddExternalForce3D<T, Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
        if (model != freeSurface) {
            modified[10] = modif::nothing;          // Fluid 2.
            modified[11] = modif::nothing;          // rhoBar2.
            modified[12] = modif::staticVariables;  // j2.
            modified[13] = modif::nothing;          // mass2.
        }
    }

private:
    T rhoDefault;
    TwoPhaseModel model;
};

template <typename T, template <typename U> class Descriptor>
struct TwoPhaseFields3D {
    static const int envelopeWidth;
    static const int smallEnvelopeWidth;
    static const int envelopeWidthForImmersedWalls;

    TwoPhaseFields3D(
        SparseBlockStructure3D const &blockStructure, Dynamics<T, Descriptor> *dynamics_,
        Dynamics<T, Descriptor> *dynamics2_, T rhoDefault_, T densityRatio_, T surfaceTension_,
        T contactAngle_, Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force2_, TwoPhaseModel model_,
        bool useImmersedWalls = false)
        // model=true --> model adapted for air bubbles, which are prevented from collapsing.
        // model=false --> model adapted for water droplets, where shearing stresses are properly
        // treated.
        :
        dynamics(dynamics_),
        dynamics2(dynamics2_),
        rhoDefault(rhoDefault_),
        densityRatio(densityRatio_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        force(force_),
        force2(force2_),
        model(model_),
        lattice(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                smallEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(), dynamics->clone()),
        lattice2(0),
        helperLists(lattice),
        mass(lattice),
        mass2(0),
        flag(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                useImmersedWalls ? envelopeWidthForImmersedWalls : envelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<int>()),
        volumeFraction((MultiBlock3D &)flag),
        curvature(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(), envelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()),
        outsideDensity((MultiBlock3D &)curvature),
        rhoBar(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()),
        rhoBar2(0),
        j(MultiBlockManagement3D(
              blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
              useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth),
          defaultMultiBlockPolicy3D().getBlockCommunicator(),
          defaultMultiBlockPolicy3D().getCombinedStatistics(),
          defaultMultiBlockPolicy3D().getMultiTensorAccess<T, 3>()),
        j2(0),
        normal((MultiBlock3D &)curvature)
    {
        initialization(blockStructure, useImmersedWalls);
    }
    // Constructor for free-surface case
    TwoPhaseFields3D(
        SparseBlockStructure3D const &blockStructure, Dynamics<T, Descriptor> *dynamics_,
        T rhoDefault_, T surfaceTension_, T contactAngle_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_, bool useImmersedWalls = false)
        // model=true --> model adapted for air bubbles, which are prevented from collapsing.
        // model=false --> model adapted for water droplets, where shearing stresses are properly
        // treated.
        :
        dynamics(dynamics_->clone()),
        dynamics2(dynamics_),
        rhoDefault(rhoDefault_),
        densityRatio(T()),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        force(force_),
        force2(force_),
        model(freeSurface),
        lattice(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                smallEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(), dynamics->clone()),
        lattice2(0),
        helperLists(lattice),
        mass(lattice),
        mass2(0),
        flag(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                useImmersedWalls ? envelopeWidthForImmersedWalls : envelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<int>()),
        volumeFraction((MultiBlock3D &)flag),
        curvature(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(), envelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()),
        outsideDensity((MultiBlock3D &)curvature),
        rhoBar(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()),
        rhoBar2(0),
        j(MultiBlockManagement3D(
              blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
              useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth),
          defaultMultiBlockPolicy3D().getBlockCommunicator(),
          defaultMultiBlockPolicy3D().getCombinedStatistics(),
          defaultMultiBlockPolicy3D().getMultiTensorAccess<T, 3>()),
        j2(0),
        normal((MultiBlock3D &)curvature)
    {
        initialization(blockStructure, useImmersedWalls);
    }
    void initialization(SparseBlockStructure3D const &blockStructure, bool useImmersedWalls)
    {
        // TwoPhaseFields3D does not work with incompressible dynamics at the moment.
        PLB_ASSERT(!dynamics->velIsJ());
        PLB_ASSERT(!dynamics2->velIsJ());

        useSurfaceTension = !util::isZero(surfaceTension);

        if (model != freeSurface) {
            lattice2 = new MultiBlockLattice3D<T, Descriptor>(
                MultiBlockManagement3D(
                    blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                    smallEnvelopeWidth),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(),
                dynamics2->clone());
            mass2 = new MultiScalarField3D<T>(lattice);
            rhoBar2 = new MultiScalarField3D<T>(
                MultiBlockManagement3D(
                    blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                    useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiScalarAccess<T>());
            j2 = new MultiTensorField3D<T, 3>(
                MultiBlockManagement3D(
                    blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                    useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiTensorAccess<T, 3>());
            lattice2->periodicity().toggleAll(true);
            mass2->periodicity().toggleAll(true);
            rhoBar2->periodicity().toggleAll(true);
            j2->periodicity().toggleAll(true);
        }

        twoPhaseArgs = aggregateTwoPhaseParams(
            lattice, lattice2, rhoBar, rhoBar2, j, j2, mass, mass2, volumeFraction, flag, normal,
            helperLists, curvature, outsideDensity);

        twoPhaseInitializeInterfaceLists3D<T, Descriptor>(helperLists);
        lattice.periodicity().toggleAll(true);
        mass.periodicity().toggleAll(true);
        flag.periodicity().toggleAll(true);
        volumeFraction.periodicity().toggleAll(true);
        curvature.periodicity().toggleAll(true);
        outsideDensity.periodicity().toggleAll(true);
        rhoBar.periodicity().toggleAll(true);
        j.periodicity().toggleAll(true);
        normal.periodicity().toggleAll(true);
        // setToConstant(flag, flag.getBoundingBox(), (int)freeSurfaceFlag::empty);
        // setToConstant(outsideDensity, outsideDensity.getBoundingBox(), rhoDefault);
        rhoBarJparam.push_back(&lattice);
        rhoBarJparam.push_back(&rhoBar);
        rhoBarJparam.push_back(&j);

        if (model != freeSurface) {
            rhoBarJparam2.push_back(lattice2);
            rhoBarJparam2.push_back(rhoBar2);
            rhoBarJparam2.push_back(j2);
        }

        lattice.internalStatSubscription().subscribeSum();     // Total mass.
        lattice.internalStatSubscription().subscribeSum();     // Lost mass.
        lattice.internalStatSubscription().subscribeIntSum();  // Num interface cells.
        if (model != freeSurface) {
            lattice2->internalStatSubscription().subscribeSum();  // Total mass 2.
            lattice2->internalStatSubscription().subscribeSum();  // Lost mass 2.
        }

        freeSurfaceDataProcessors();
    }

    TwoPhaseFields3D(TwoPhaseFields3D<T, Descriptor> const &rhs) :
        dynamics(rhs.dynamics->clone()),
        dynamics2(rhs.dynamics2->clone()),
        rhoDefault(rhs.rhoDefault),
        densityRatio(rhs.densityRatio),
        surfaceTension(rhs.surfaceTension),
        contactAngle(rhs.contactAngle),
        useSurfaceTension(rhs.useSurfaceTension),
        force(rhs.force),
        force2(rhs.force2),
        model(rhs.model),
        lattice(rhs.lattice),
        lattice2(rhs.lattice2 ? new MultiBlockLattice3D<T, Descriptor>(*rhs.lattice2) : 0),
        helperLists(rhs.helperLists),
        mass(rhs.mass),
        mass2(rhs.mass2 ? new MultiScalarField3D<T>(*rhs.mass2) : 0),
        flag(rhs.flag),
        volumeFraction(rhs.volumeFraction),
        curvature(rhs.curvature),
        outsideDensity(rhs.outsideDensity),
        rhoBar(rhs.rhoBar),
        rhoBar2(rhs.rhoBar2 ? new MultiScalarField3D<T>(*rhs.rhoBar2) : 0),
        j(rhs.j),
        j2(rhs.j2 ? new MultiTensorField3D<T, 3>(*rhs.j2) : 0),
        normal(rhs.normal),
        rhoBarJparam(),
        rhoBarJparam2(),
        twoPhaseArgs()
    {
        twoPhaseArgs = aggregateTwoPhaseParams(
            lattice, lattice2, rhoBar, rhoBar2, j, j2, mass, mass2, volumeFraction, flag, normal,
            helperLists, curvature, outsideDensity);

        rhoBarJparam.push_back(&lattice);
        rhoBarJparam.push_back(&rhoBar);
        rhoBarJparam.push_back(&j);

        if (model != freeSurface) {
            rhoBarJparam2.push_back(&lattice2);
            rhoBarJparam2.push_back(&rhoBar2);
            rhoBarJparam2.push_back(&j2);
        }
    }

    void swap(TwoPhaseFields3D<T, Descriptor> &rhs)
    {
        std::swap(dynamics, rhs.dynamics);
        std::swap(dynamics2, rhs.dynamics2);
        std::swap(rhoDefault, rhs.rhoDefault);
        std::swap(densityRatio, rhs.densityRatio);
        std::swap(surfaceTension, rhs.surfaceTension);
        std::swap(contactAngle, rhs.contactAngle);
        std::swap(useSurfaceTension, rhs.useSurfaceTension);
        std::swap(force, rhs.force);
        std::swap(force2, rhs.force2);
        std::swap(model, rhs.model);
        std::swap(lattice, rhs.lattice);
        std::swap(lattice2, rhs.lattice2);
        std::swap(helperLists, rhs.helperLists);
        std::swap(mass, rhs.mass);
        std::swap(mass2, rhs.mass2);
        std::swap(flag, rhs.flag);
        std::swap(volumeFraction, rhs.volumeFraction);
        std::swap(curvature, rhs.curvature);
        std::swap(outsideDensity, rhs.outsideDensity);
        std::swap(rhoBar, rhs.rhoBar);
        std::swap(rhoBar2, rhs.rhoBar2);
        std::swap(j2, rhs.j2);
        std::swap(normal, rhs.normal);
        std::swap(rhoBarJparam, rhs.rhoBarJparam);
        std::swap(rhoBarJparam2, rhs.rhoBarJparam2);
        std::swap(twoPhaseArgs, rhs.twoPhaseArgs);
    }

    TwoPhaseFields3D<T, Descriptor> &operator=(TwoPhaseFields3D<T, Descriptor> const &rhs)
    {
        TwoPhaseFields3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }

    TwoPhaseFields3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseFields3D<T, Descriptor>(*this);
    }

    ~TwoPhaseFields3D()
    {
        delete dynamics;
        delete dynamics2;
        delete lattice2;
        delete rhoBar2;
        delete j2;
        delete mass2;
    }

    void periodicityToggle(plint direction, bool periodic)
    {
        PLB_ASSERT(direction == 0 || direction == 1 || direction == 2);

        lattice.periodicity().toggle(direction, periodic);
        mass.periodicity().toggle(direction, periodic);
        flag.periodicity().toggle(direction, periodic);
        volumeFraction.periodicity().toggle(direction, periodic);
        curvature.periodicity().toggle(direction, periodic);
        outsideDensity.periodicity().toggle(direction, periodic);
        rhoBar.periodicity().toggle(direction, periodic);
        j.periodicity().toggle(direction, periodic);
        normal.periodicity().toggle(direction, periodic);
        if (model != freeSurface) {
            lattice2->periodicity().toggle(direction, periodic);
            mass2->periodicity().toggle(direction, periodic);
            rhoBar2->periodicity().toggle(direction, periodic);
            j2->periodicity().toggle(direction, periodic);
        }
    }

    void periodicityToggleAll(bool periodic)
    {
        lattice.periodicity().toggleAll(periodic);
        mass.periodicity().toggleAll(periodic);
        flag.periodicity().toggleAll(periodic);
        volumeFraction.periodicity().toggleAll(periodic);
        curvature.periodicity().toggleAll(periodic);
        outsideDensity.periodicity().toggleAll(periodic);
        rhoBar.periodicity().toggleAll(periodic);
        j.periodicity().toggleAll(periodic);
        normal.periodicity().toggleAll(periodic);
        if (model != freeSurface) {
            lattice2->periodicity().toggleAll(periodic);
            mass2->periodicity().toggleAll(periodic);
            rhoBar2->periodicity().toggleAll(periodic);
            j2->periodicity().toggleAll(periodic);
        }
    }

    void defaultInitialize()
    {
        applyProcessingFunctional(
            new DefaultInitializeTwoPhase3D<T, Descriptor>(
                dynamics->clone(), dynamics2->clone(), force, force2, rhoDefault,
                model == freeSurface),
            lattice.getBoundingBox(), twoPhaseArgs);
    }

    void partiallyDefaultInitialize()
    {
        applyProcessingFunctional(
            new PartiallyDefaultInitializeTwoPhase3D<T, Descriptor>(
                dynamics->clone(), dynamics2->clone(), force, force2, rhoDefault,
                model == freeSurface),
            lattice.getBoundingBox(), twoPhaseArgs);
    }
    std::unique_ptr<MultiScalarField3D<T> > computePressure(
        Box3D domain, bool computeFluid1, bool computeFluid2)
    {
        std::unique_ptr<MultiScalarField3D<T> > pressure =
            generateMultiScalarField<T>(lattice, domain);
        std::vector<MultiBlock3D *> args(twoPhaseArgs);
        args.push_back(pressure.get());
        applyProcessingFunctional(
            new TwoPhaseComputePressure3D<T, Descriptor>(
                densityRatio, rhoDefault, model, computeFluid1, computeFluid2),
            domain, args);
        return pressure;
    }
    std::unique_ptr<MultiScalarField3D<T> > computePressure(Box3D domain)
    {
        return computePressure(domain, true, true);
    }
    std::unique_ptr<MultiScalarField3D<T> > computePressure()
    {
        return computePressure(lattice.getBoundingBox(), true, true);
    }
    void computePressureAverage(Box3D domain, T &p1, T &p2)
    {
        TwoPhaseAveragePressure3D<T, Descriptor> functional(densityRatio, rhoDefault, model);
        applyProcessingFunctional(functional, domain, twoPhaseArgs);
        p1 = functional.getAveragePressure();
        p2 = functional.getAveragePressure2();
    }
    void computeVelocityAverage(Box3D domain, Array<T, 3> &v1, Array<T, 3> &v2)
    {
        TwoPhaseAverageVelocity3D<T, Descriptor> functional(model);
        applyProcessingFunctional(functional, domain, twoPhaseArgs);
        v1 = functional.getAverageVelocity();
        v2 = functional.getAverageVelocity2();
    }
    std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocity(
        Box3D domain, bool computeFluid1, bool computeFluid2)
    {
        std::unique_ptr<MultiTensorField3D<T, 3> > velocity =
            generateMultiTensorField<T, 3>(lattice, domain);
        std::vector<MultiBlock3D *> args(twoPhaseArgs);
        args.push_back(velocity.get());
        applyProcessingFunctional(
            new TwoPhaseComputeVelocity3D<T, Descriptor>(
                densityRatio, computeFluid1, computeFluid2, model == freeSurface),
            domain, args);
        return velocity;
    }
    std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocity(Box3D domain)
    {
        return computeVelocity(domain, true, true);
    }
    std::unique_ptr<MultiTensorField3D<T, 3> > computeVelocity()
    {
        return computeVelocity(lattice.getBoundingBox(), true, true);
    }

    void freeSurfaceDataProcessors()
    {
        plint pl;  // Processor level.

        /***** Initial level ******/
        pl = 0;

        integrateProcessingFunctional(
            new ExternalRhoJcollideAndStream3D<T, Descriptor>, lattice.getBoundingBox(),
            rhoBarJparam, pl);

        // twophase
        if (model != freeSurface) {
            integrateProcessingFunctional(
                new ExternalRhoJcollideAndStream3D<T, Descriptor>, lattice.getBoundingBox(),
                rhoBarJparam2, pl);
        }

        integrateProcessingFunctional(
            new FreeSurfaceComputeNormals3D<T, Descriptor>, lattice.getBoundingBox(), twoPhaseArgs,
            pl);

        /***** New level ******/
        pl++;

        if (useSurfaceTension) {
            integrateProcessingFunctional(
                new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngle),
                lattice.getBoundingBox(), twoPhaseArgs, pl);

            // To change to the curvature calculation with height functions, uncomment the next data
            // processor and comment out the two previous ones. If only the next data processor is
            // used and there is no surface tension, the normals are not computed at all. Be careful
            // if you intent to use the normals and do not have the surface tension algorithm
            // enabled.
            // integrateProcessingFunctional (
            //        new FreeSurfaceGeometry3D<T,Descriptor>(contactAngle),
            //        lattice.getBoundingBox(), freeSurfaceArgs, pl );
        }

        integrateProcessingFunctional(
            new TwoPhaseMassChange3D<T, Descriptor>(model == freeSurface), lattice.getBoundingBox(),
            twoPhaseArgs, pl);

        integrateProcessingFunctional(
            new TwoPhaseCompletion3D<T, Descriptor>(model == freeSurface), lattice.getBoundingBox(),
            twoPhaseArgs, pl);

        integrateProcessingFunctional(
            new TwoPhaseMacroscopic3D<T, Descriptor>(
                rhoDefault, densityRatio, surfaceTension, model),
            lattice.getBoundingBox(), twoPhaseArgs, pl);

        /***** New level ******/
        // pl++;

        // integrateProcessingFunctional (
        //         new TwoPhaseInterfaceFilter<T,Descriptor>(rhoDefault, densityRatio, model),
        //         lattice.getBoundingBox(), twoPhaseArgs, pl );

        /***** New level ******/
        // pl++;

        // integrateProcessingFunctional (
        //         new TwoPhaseInterfaceFilter<T,Descriptor>(rhoDefault, densityRatio, model),
        //         lattice.getBoundingBox(), twoPhaseArgs, pl );

        if (useSurfaceTension) {
            // TwoPhaseFields3D does not work with incompressible dynamics at the moment.
            bool incompressibleModel = false;
            integrateProcessingFunctional(
                new FreeSurfaceAddSurfaceTension3D<T, Descriptor>(
                    surfaceTension, incompressibleModel),
                lattice.getBoundingBox(), twoPhaseArgs, pl);
        }

        integrateProcessingFunctional(
            new TwoPhaseStabilize3D<T, Descriptor>(model), lattice.getBoundingBox(), twoPhaseArgs,
            pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new TwoPhaseComputeInterfaceLists3D<T, Descriptor>(
                model, rhoDefault, densityRatio, surfaceTension),
            lattice.getBoundingBox(), twoPhaseArgs, pl);

        // interface->fluid   --     interface->empty
        // interface->empty   --     interface->fluid
        // fluid->interface   --     empty->interface
        integrateProcessingFunctional(
            new TwoPhaseIniInterfaceToAnyNodes3D<T, Descriptor>(
                rhoDefault, dynamics->clone(), dynamics2->clone(), force, force2, densityRatio,
                surfaceTension, model),
            lattice.getBoundingBox(), twoPhaseArgs, pl);

        // empty->interface  --      fluid->interface
        integrateProcessingFunctional(
            new TwoPhaseIniEmptyToInterfaceNodes3D<T, Descriptor>(
                dynamics->clone(), force, densityRatio, model),
            lattice.getBoundingBox(), twoPhaseArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new TwoPhaseRemoveFalseInterfaceCells3D<T, Descriptor>(rhoDefault, model),
            lattice.getBoundingBox(), twoPhaseArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new TwoPhaseEqualMassExcessReDistribution3D<T, Descriptor>(model),
            lattice.getBoundingBox(), twoPhaseArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceComputeStatistics3D<T, Descriptor>, lattice.getBoundingBox(),
            twoPhaseArgs, pl);

        bool useForce = !util::isZero(norm(force));
        bool useForce2 = !util::isZero(norm(force2));
        if (useForce || useForce2) {
            integrateProcessingFunctional(
                new TwoPhaseAddExternalForce3D<T, Descriptor>(rhoDefault, model),
                lattice.getBoundingBox(), twoPhaseArgs, pl);
        }

        // integrateProcessingFunctional (
        //     new VerifyTwoPhase<T,Descriptor>(),
        //     lattice.getBoundingBox(), twoPhaseArgs, pl);
    }

    void appendBlocksToCheckpointVector(std::vector<MultiBlock3D *> &checkpointBlocks)
    {
        checkpointBlocks.push_back(&lattice);
        checkpointBlocks.push_back(&mass);
        checkpointBlocks.push_back(&flag);
        checkpointBlocks.push_back(&volumeFraction);
        checkpointBlocks.push_back(&outsideDensity);
        checkpointBlocks.push_back(&rhoBar);
        checkpointBlocks.push_back(&j);
        if (model != freeSurface) {
            checkpointBlocks.push_back(lattice2);
            checkpointBlocks.push_back(mass2);
            checkpointBlocks.push_back(rhoBar2);
            checkpointBlocks.push_back(j2);
        }
    }

    Dynamics<T, Descriptor> *dynamics, *dynamics2;
    T rhoDefault;
    T densityRatio;
    T surfaceTension;
    T contactAngle;
    bool useSurfaceTension;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force, force2;
    TwoPhaseModel model;
    MultiBlockLattice3D<T, Descriptor> lattice;
    MultiBlockLattice3D<T, Descriptor> *lattice2;
    MultiContainerBlock3D helperLists;
    MultiScalarField3D<T> mass;
    MultiScalarField3D<T> *mass2;
    MultiScalarField3D<int> flag;
    MultiScalarField3D<T> volumeFraction;
    MultiScalarField3D<T> curvature;
    MultiScalarField3D<T> outsideDensity;
    MultiScalarField3D<T> rhoBar;
    MultiScalarField3D<T> *rhoBar2;
    MultiTensorField3D<T, 3> j;
    MultiTensorField3D<T, 3> *j2;
    MultiTensorField3D<T, 3> normal;
    std::vector<MultiBlock3D *> rhoBarJparam, rhoBarJparam2;
    std::vector<MultiBlock3D *> twoPhaseArgs;
};

template <typename T, template <typename U> class Descriptor>
const int TwoPhaseFields3D<T, Descriptor>::envelopeWidth =
    3;  // Necessary when we use height functions to compute the curvature,
        // or when double smoothing is used at the data processor that
        // computes the normals from the volume fraction.
// template<typename T, template<typename U> class Descriptor>
// const int TwoPhaseFields3D<T,Descriptor>::envelopeWidth = 4; // Necessary when we use height
// functions to compute the curvature and
//  use the old contact angle algorithm.
template <typename T, template <typename U> class Descriptor>
const int TwoPhaseFields3D<T, Descriptor>::smallEnvelopeWidth = 1;

template <typename T, template <typename U> class Descriptor>
const int TwoPhaseFields3D<T, Descriptor>::envelopeWidthForImmersedWalls = 4;

template <typename T, template <typename U> class Descriptor>
class TwoPhaseOutletMaximumVolumeFraction3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseOutletMaximumVolumeFraction3D(T volumeFraction_, TwoPhaseModel model_) :
        volumeFraction(volumeFraction_), model(model_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TwoPhaseOutletMaximumVolumeFraction3D<T, Descriptor> *clone() const
    {
        return new TwoPhaseOutletMaximumVolumeFraction3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        PLB_ASSERT(modified.size() >= 14);
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables;   // Fluid.
        modified[1] = modif::staticVariables;   // rhoBar.
        modified[2] = modif::nothing;           // j.
        modified[3] = modif::staticVariables;   // Mass.
        modified[4] = modif::nothing;           // Volume fraction.
        modified[5] = modif::nothing;           // Flag-status.
        modified[6] = modif::nothing;           // Normal.
        modified[7] = modif::nothing;           // Interface lists.
        modified[8] = modif::nothing;           // Curvature.
        modified[9] = modif::nothing;           // Outside density.
        modified[10] = modif::staticVariables;  // Fluid 2.
        modified[11] = modif::staticVariables;  // rhoBar2.
        modified[12] = modif::staticVariables;  // j2.
        modified[13] = modif::staticVariables;  // mass2.
    }

private:
    T volumeFraction;
    TwoPhaseModel model;
};

}  // namespace plb

#endif  // TWO_PHASE_MODEL_3D_H
