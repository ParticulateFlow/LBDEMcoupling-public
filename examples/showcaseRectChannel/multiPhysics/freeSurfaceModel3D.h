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

#ifndef FREE_SURFACE_MODEL_3D_H
#define FREE_SURFACE_MODEL_3D_H

#include <algorithm>
#include <map>
#include <set>

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "basicDynamics/dynamicsProcessor3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiBlock/coupling3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/group3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "multiPhysics/freeSurfaceInitializer3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"
#include "offLattice/immersedWalls3D.h"

#ifdef PLB_MPI_PARALLEL
// DISABLE_WARNING_PUSH
// DISABLE_WARNING_CAST_FUNCTION_TYPE
#include <mpi.h>
// DISABLE_WARNING_POP
#endif

namespace plb {

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceComputeNormals3D : public BoxProcessingFunctional3D {
public:
    virtual FreeSurfaceComputeNormals3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceComputeNormals3D<T, Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::staticVariables;  // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }
};

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceGeometry3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceGeometry3D(T contactAngle_) : contactAngle(contactAngle_)
    {
        // The contact angle must take values between 0 and 180 degrees. If it is negative,
        // this means that contact angle effects will not be modeled.
        PLB_ASSERT(util::lessEqual(contactAngle, (T)180));

        if (util::lessThan(contactAngle, (T)0)) {
            useContactAngle = 0;
        } else {
            useContactAngle = 1;
        }
    }
    virtual FreeSurfaceGeometry3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceGeometry3D<T, Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::staticVariables;  // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::staticVariables;  // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    ScalarField3D<int> *getInterfaceFlags(
        Box3D domain, FreeSurfaceProcessorParam3D<T, Descriptor> &param);
    void computeHeights3D(
        FreeSurfaceProcessorParam3D<T, Descriptor> &param, int integrationDirection, plint iX,
        plint iY, plint iZ, T h[3][3]);
    void computeHeights2D(
        FreeSurfaceProcessorParam3D<T, Descriptor> &param, Array<int, 3> &wallTangent0,
        Array<int, 3> &wallTangent1, int integrationDirection, plint iX, plint iY, plint iZ,
        T h[3]);

private:
    enum { unTagged = 0, notInterface = 1, regular = 2, contactLine = 4, adjacent = 10 };

private:
    T contactAngle;
    int useContactAngle;
};

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceComputeCurvature3D : public BoxProcessingFunctional3D {
public:
    typedef T (*ContactAngleFunction)(T x, T y, T z);  // Returns the contact angle in degrees.
public:
    FreeSurfaceComputeCurvature3D(T contactAngle_) :
        contactAngle(contactAngle_), contactAngleFunction(0)
    {
        // The contact angle must take values between 0 and 180 degrees. If it is negative,
        // this means that contact angle effects will not be modeled.
        PLB_ASSERT(util::lessEqual(contactAngle, (T)180));

        if (util::lessThan(contactAngle, (T)0)) {
            useContactAngle = 0;
        } else {
            useContactAngle = 1;
        }

        if (useContactAngle) {
            T pi = std::acos((T)-1);
            contactAngle *= pi / (T)180;
        }
    }
    FreeSurfaceComputeCurvature3D(ContactAngleFunction contactAngleFunction_) :
        contactAngle(-1.0), contactAngleFunction(contactAngleFunction_)
    {
        // The contact angle must take values between 0 and 180 degrees.
        // If the function pointer is 0, this means that contact angle effects will not be modeled.
        if (contactAngleFunction == 0) {
            useContactAngle = 0;
        } else {
            useContactAngle = 1;
        }
    }
    virtual FreeSurfaceComputeCurvature3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceComputeCurvature3D<T, Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid.
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::nothing;          // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::staticVariables;  // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    T contactAngle;
    ContactAngleFunction contactAngleFunction;
    int useContactAngle;
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
class FreeSurfaceMassChange3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceMassChange3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceMassChange3D<T, Descriptor>(*this);
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
    }
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
class FreeSurfaceCompletion3D : public BoxProcessingFunctional3D {
public:
    virtual FreeSurfaceCompletion3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceCompletion3D<T, Descriptor>(*this);
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
    }
};

/// Compute and store mass, volume-fraction and macroscopic variables.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Mass:          needed in bulk
 *   - Populations:   needed in bulk
 * Output:
 *   - mass, volume-fraction, density, momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceMacroscopic3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceMacroscopic3D(bool incompressibleModel_) : incompressibleModel(incompressibleModel_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceMacroscopic3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceMacroscopic3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid. Should be: staticVariables (maybe not...).
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass. Should be: staticVariables.
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status. TODO: used to be staticVariables...
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    bool incompressibleModel;
};

/// The same as FreeSurfaceMacroscopic3D, but no global statistics from the previous time
/// step are used for lost mass redistribution.
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceMacroscopicWithoutLostMassReDistribution3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceMacroscopicWithoutLostMassReDistribution3D(bool incompressibleModel_) :
        incompressibleModel(incompressibleModel_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceMacroscopicWithoutLostMassReDistribution3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceMacroscopicWithoutLostMassReDistribution3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid. Should be: staticVariables (maybe not...).
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::nothing;          // Mass.
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status. TODO: used to be staticVariables...
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    bool incompressibleModel;
};

/// Add the surface tension contribution.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Mass:          needed in bulk
 *   - Populations:   needed in bulk
 * Output:
 *   - volume-fraction, density, momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceAddSurfaceTension3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceAddSurfaceTension3D(T surfaceTension_, bool incompressibleModel_) :
        surfaceTension(surfaceTension_), incompressibleModel(incompressibleModel_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceAddSurfaceTension3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceAddSurfaceTension3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid. Should be: staticVariables (maybe not...).
        modified[1] = modif::staticVariables;  // rhoBar.
        if (incompressibleModel) {
            modified[2] = modif::nothing;  // j.
        } else {
            modified[2] = modif::staticVariables;  // j.
        }
        modified[3] = modif::nothing;          // Mass. TODO: used to be staticVariables...
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status. TODO: used to be staticVariables...
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    T surfaceTension;
    bool incompressibleModel;
};

/// Add the surface tension contribution by using a Weber number model.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Mass:          needed in bulk
 *   - Populations:   needed in bulk
 * Output:
 *   - volume-fraction, density, momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceAddSurfaceTensionWeberModel3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceAddSurfaceTensionWeberModel3D(
        T surfaceTension_, T rhoDefault_, T characteristicLength_, T criticalWe_, T criticalKn_,
        bool incompressibleModel_) :
        surfaceTension(surfaceTension_),
        rhoDefault(rhoDefault_),
        characteristicLength(characteristicLength_),
        criticalWe(criticalWe_),
        criticalKn(criticalKn_),
        incompressibleModel(incompressibleModel_)
    {
        PLB_ASSERT(surfaceTension > (T)0);
        PLB_ASSERT(rhoDefault > (T)0);
        PLB_ASSERT(characteristicLength > (T)0);
        PLB_ASSERT(criticalWe > (T)0);
        PLB_ASSERT(criticalKn > (T)0);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceAddSurfaceTensionWeberModel3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceAddSurfaceTensionWeberModel3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid. Should be: staticVariables (maybe not...).
        modified[1] = modif::staticVariables;  // rhoBar.
        if (incompressibleModel) {
            modified[2] = modif::nothing;  // j.
        } else {
            modified[2] = modif::staticVariables;  // j.
        }
        modified[3] = modif::nothing;          // Mass. TODO: used to be staticVariables...
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status. TODO: used to be staticVariables...
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    T surfaceTension, rhoDefault, characteristicLength, criticalWe, criticalKn;
    bool incompressibleModel;
};

/// Add the surface tension contribution. The surface tension coefficient is read from a scalar
/// field.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Mass:          needed in bulk
 *   - Populations:   needed in bulk
 * Output:
 *   - volume-fraction, density, momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceAddSurfaceTensionFromScalarField3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceAddSurfaceTensionFromScalarField3D(bool incompressibleModel_) :
        incompressibleModel(incompressibleModel_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceAddSurfaceTensionFromScalarField3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceAddSurfaceTensionFromScalarField3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid. Should be: staticVariables (maybe not...).
        modified[1] = modif::staticVariables;  // rhoBar.
        if (incompressibleModel) {
            modified[2] = modif::nothing;  // j.
        } else {
            modified[2] = modif::staticVariables;  // j.
        }
        modified[3] = modif::nothing;          // Mass. TODO: used to be staticVariables...
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status. TODO: used to be staticVariables...
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface-lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.

        modified[10] = modif::nothing;  // Surface tension coefficient.
    }

private:
    bool incompressibleModel;
};

/// Stabilization scheme on the post-collide populations.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Populations:   needed in bulk
 *   - Momentum:      needed in bulk
 *   - Density:       needed in bulk
 * Output:
 *   - Populations.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceStabilize3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceStabilize3D() { }
    virtual FreeSurfaceStabilize3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceStabilize3D<T, Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;  // Fluid (should be staticVariables).
        modified[1] = modif::nothing;  // rhoBar.
        modified[2] = modif::nothing;  // j.
        modified[3] = modif::nothing;  // Mass.
        modified[4] = modif::nothing;  // Volume fraction.
        modified[5] = modif::nothing;  // Flag-status.
        modified[6] = modif::nothing;  // Normal.
        modified[7] = modif::nothing;  // Interface lists.
        modified[8] = modif::nothing;  // Curvature.
        modified[9] = modif::nothing;  // Outside density.
    }
};

/// Based on the current flag status, decide, upon the value of mass fraction, which nodes shall
///   switch state.
/** Input:
 *   - Volume fraction: needed in bulk+2
 *   - Flag-status:   needed in bulk+2
 * Output:
 *   - interface-to-fluid list: defined in bulk+2
 *   - interface-to-empty list: defined in bulk+1
 *   - empty-to-interface list: defined in bulk+1
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceComputeInterfaceLists3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceComputeInterfaceLists3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid (not used in this processor).
        modified[1] = modif::nothing;          // rhoBar.
        modified[2] = modif::nothing;          // j.
        modified[3] = modif::nothing;          // Mass (not used in this processor).
        modified[4] = modif::nothing;          // Volume fraction, read-only.
        modified[5] = modif::nothing;          // Flag-status, read-only.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::staticVariables;  // Interface-lists; all lists are created here.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    static T kappa;  // Safety threshold for state-change, to prevent back-and-forth oscillations.
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
class FreeSurfaceIniInterfaceToAnyNodes3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceIniInterfaceToAnyNodes3D(
        T rhoDefault_, Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate_) :
        rhoDefault(rhoDefault_), emptyNodeDynamicsTemplate(emptyNodeDynamicsTemplate_)
    { }
    FreeSurfaceIniInterfaceToAnyNodes3D(T rhoDefault_) : rhoDefault(rhoDefault_)
    {
        emptyNodeDynamicsTemplate = new NoDynamics<T, Descriptor>(rhoDefault);
    }
    FreeSurfaceIniInterfaceToAnyNodes3D(
        FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor> const &rhs) :
        rhoDefault(rhs.rhoDefault),
        emptyNodeDynamicsTemplate(rhs.emptyNodeDynamicsTemplate->clone())
    { }
    FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor> *operator=(
        FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor> const &rhs)
    {
        FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor> &rhs)
    {
        std::swap(rhoDefault, rhs.rhoDefault);
        std::swap(emptyNodeDynamicsTemplate, rhs.emptyNodeDynamicsTemplate);
    }
    ~FreeSurfaceIniInterfaceToAnyNodes3D()
    {
        delete emptyNodeDynamicsTemplate;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(*this);
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
    }

private:
    T rhoDefault;
    Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate;
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
class FreeSurfaceIniEmptyToInterfaceNodes3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceIniEmptyToInterfaceNodes3D(
        Dynamics<T, Descriptor> *dynamicsTemplate_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_) :
        dynamicsTemplate(dynamicsTemplate_), force(force_)
    { }
    FreeSurfaceIniEmptyToInterfaceNodes3D(
        FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor> const &rhs) :
        dynamicsTemplate(rhs.dynamicsTemplate->clone()), force(rhs.force)
    { }
    FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor> *operator=(
        FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor> const &rhs)
    {
        FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor> &rhs)
    {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(force, rhs.force);
    }
    ~FreeSurfaceIniEmptyToInterfaceNodes3D()
    {
        delete dynamicsTemplate;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(*this);
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
    }

private:
    Dynamics<T, Descriptor> *dynamicsTemplate;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce>
        force;  // Body force, for initialization of the new interface cell.
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
class FreeSurfaceRemoveFalseInterfaceCells3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceRemoveFalseInterfaceCells3D(
        T rhoDefault_, Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate_) :
        rhoDefault(rhoDefault_), emptyNodeDynamicsTemplate(emptyNodeDynamicsTemplate_)
    { }
    FreeSurfaceRemoveFalseInterfaceCells3D(T rhoDefault_) : rhoDefault(rhoDefault_)
    {
        emptyNodeDynamicsTemplate = new NoDynamics<T, Descriptor>(rhoDefault);
    }
    FreeSurfaceRemoveFalseInterfaceCells3D(
        FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor> const &rhs) :
        rhoDefault(rhs.rhoDefault),
        emptyNodeDynamicsTemplate(rhs.emptyNodeDynamicsTemplate->clone())
    { }
    FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor> *operator=(
        FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor> const &rhs)
    {
        FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor> &rhs)
    {
        std::swap(rhoDefault, rhs.rhoDefault);
        std::swap(emptyNodeDynamicsTemplate, rhs.emptyNodeDynamicsTemplate);
    }
    ~FreeSurfaceRemoveFalseInterfaceCells3D()
    {
        delete emptyNodeDynamicsTemplate;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(*this);
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
    }

private:
    T rhoDefault;
    Dynamics<T, Descriptor> *emptyNodeDynamicsTemplate;
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
class FreeSurfaceEqualMassExcessReDistribution3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceEqualMassExcessReDistribution3D(*this);
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
    }
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
class FreeSurfaceEqualMassExcessReDistributionAndComputationOfLostMass3D :
    public BoxProcessingFunctional3D {
public:
    FreeSurfaceEqualMassExcessReDistributionAndComputationOfLostMass3D(
        FreeSurfaceReductionData &reductionData_) :
        reductionData(reductionData_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceEqualMassExcessReDistributionAndComputationOfLostMass3D<T, Descriptor>
        *clone() const
    {
        return new FreeSurfaceEqualMassExcessReDistributionAndComputationOfLostMass3D(*this);
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
    }

private:
    FreeSurfaceReductionData &reductionData;
};

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceComputeStatistics3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceComputeStatistics3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceComputeStatistics3D(*this);
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
    }
};

template <typename T>
class FreeSurfaceComputeReductionsPerProcess3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceComputeReductionsPerProcess3D(FreeSurfaceReductionData &reductionData_) :
        reductionData(reductionData_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;  // Flag.
        modified[1] = modif::nothing;  // Mass.
    }
    virtual FreeSurfaceComputeReductionsPerProcess3D<T> *clone() const
    {
        return new FreeSurfaceComputeReductionsPerProcess3D<T>(*this);
    }

private:
    FreeSurfaceReductionData &reductionData;
};

#ifdef PLB_MPI_PARALLEL
// CAUTION: This data processor will not work (will hang) if the multi-block passed does
//          not have at least one atomic block on every MPI process.
class FreeSurfaceComputeReductions3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceComputeReductions3D(
        FreeSurfaceReductionData &reductionData_, MPI_Comm &reductionCommunicator_) :
        reductionData(reductionData_), reductionCommunicator(reductionCommunicator_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    virtual FreeSurfaceComputeReductions3D *clone() const
    {
        return new FreeSurfaceComputeReductions3D(*this);
    }

private:
    FreeSurfaceReductionData &reductionData;
    MPI_Comm &reductionCommunicator;
};
#else
class FreeSurfaceComputeSerialReductions3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceComputeSerialReductions3D(FreeSurfaceReductionData &reductionData_) :
        reductionData(reductionData_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    virtual FreeSurfaceComputeSerialReductions3D *clone() const
    {
        return new FreeSurfaceComputeSerialReductions3D(*this);
    }

private:
    FreeSurfaceReductionData &reductionData;
};
#endif

class FreeSurfaceResetReductionData3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceResetReductionData3D(FreeSurfaceReductionData &reductionData_) :
        reductionData(reductionData_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> blocks);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    virtual FreeSurfaceResetReductionData3D *clone() const
    {
        return new FreeSurfaceResetReductionData3D(*this);
    }

private:
    FreeSurfaceReductionData &reductionData;
};

template <typename T, template <typename U> class Descriptor>
class InitializeInterfaceLists3D : public BoxProcessingFunctional3D {
    virtual void processGenericBlocks(
        [[maybe_unused]] Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
    {
        PLB_ASSERT(atomicBlocks.size() == 1);

        AtomicContainerBlock3D *containerInterfaceLists =
            dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[0]);
        PLB_ASSERT(containerInterfaceLists);
        InterfaceLists<T, Descriptor> *interfaceLists = new InterfaceLists<T, Descriptor>;
        containerInterfaceLists->setData(interfaceLists);
    }
    virtual InitializeInterfaceLists3D<T, Descriptor> *clone() const
    {
        return new InitializeInterfaceLists3D<T, Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        // Default-assign potential other parameters present in a multi-fluid system.
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables;
    }
};

/// Wrapper for execution of InitializeInterfaceLists3D.
template <typename T, template <typename U> class Descriptor>
void initializeInterfaceLists3D(MultiContainerBlock3D &interfaceListBlock)
{
    std::vector<MultiBlock3D *> arg;
    arg.push_back(&interfaceListBlock);
    applyProcessingFunctional(
        new InitializeInterfaceLists3D<T, Descriptor>, interfaceListBlock.getBoundingBox(), arg);
}

/// Addition of the external forces.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Momentum:      needed in bulk
 * Output:
 *   - Momentum.
 **/
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceAddExternalForce3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceAddExternalForce3D(T rhoDefault_) : rhoDefault(rhoDefault_) { }
    virtual FreeSurfaceAddExternalForce3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceAddExternalForce3D<T, Descriptor>(*this);
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
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }

private:
    T rhoDefault;
};

/// Re-distribute lost mass if the FreeSurfaceMacroscopicWithoutLostMassReDistribution3D is used
/// instead of the FreeSurfaceMacroscopic3D one.
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceLostMassReDistribution3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceLostMassReDistribution3D(FreeSurfaceReductionData &reductionData_) :
        reductionData(reductionData_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual FreeSurfaceLostMassReDistribution3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceLostMassReDistribution3D<T, Descriptor>(*this);
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
    }

private:
    FreeSurfaceReductionData &reductionData;
};

/// Repelling interface nodes kinematically from the immersed walls.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Density:       needed in bulk
 *   - Momentum:      needed in bulk
 *   - Container:     needed in bulk and envelope
 * Output:
 *   - Density.
 *   - Momentum.
 **/
// TODO: This data processor changes the momentum in the vicinity of the
//       immersed boundary. Maybe it affects strongly the measurement of
//       force and torque on the immersed surface. Needs to be checked.
template <typename T, class VelFunction>
class RepelInterfaceFromImmersedWalls3D : public BoxProcessingFunctional3D {
public:
    RepelInterfaceFromImmersedWalls3D(
        VelFunction velFunction_, T rhoDefault_, bool strongRepelling_) :
        velFunction(velFunction_), rhoDefault(rhoDefault_), strongRepelling(strongRepelling_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual RepelInterfaceFromImmersedWalls3D<T, VelFunction> *clone() const
    {
        return new RepelInterfaceFromImmersedWalls3D<T, VelFunction>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;  // RhoBar
        modified[1] = modif::staticVariables;  // J
        modified[2] = modif::nothing;          // Flag
        modified[3] = modif::nothing;          // Container Block with triangle data.
    }
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }

private:
    VelFunction velFunction;
    T rhoDefault;
    bool strongRepelling;
};

/// Protect temporarily fluid nodes close to immersed walls from turing to interface.
/** Input:
 *   - Flag-status:   needed in bulk
 *   - Container:     needed in bulk and envelope
 * Output:
 *   - Flag-status.
 **/
template <typename T>
class TemporarilyProtectImmersedWalls3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual TemporarilyProtectImmersedWalls3D<T> *clone() const
    {
        return new TemporarilyProtectImmersedWalls3D<T>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;  // Flag
        modified[1] = modif::nothing;          // Container Block with triangle data.
    }
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }
};

/// Remove the above protection.
/** Input:
 *   - Flag-status:   needed in bulk
 * Output:
 *   - Flag-status.
 **/
template <typename T>
class RemoveProtectionFromImmersedWalls3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual RemoveProtectionFromImmersedWalls3D<T> *clone() const
    {
        return new RemoveProtectionFromImmersedWalls3D<T>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::staticVariables;  // Flag
    }
    virtual BlockDomain::DomainT appliesTo() const
    {
        return BlockDomain::bulk;
    }
};

template <typename T, template <typename U> class Descriptor>
struct FreeSurfaceFields3D {
    static const int envelopeWidth;
    static const int smallEnvelopeWidth;
    static const int envelopeWidthForImmersedWalls;

    FreeSurfaceFields3D(
        SparseBlockStructure3D const &blockStructure, Dynamics<T, Descriptor> *dynamics_,
        T rhoDefault_, T surfaceTension_, T contactAngle_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_) :
        dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        force(force_),
        lattice(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                smallEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(), dynamics->clone()),
        helperLists(lattice),
        mass(lattice),
        flag(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(), envelopeWidth),
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
                smallEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()),
        j(MultiBlockManagement3D(
              blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
              smallEnvelopeWidth),
          defaultMultiBlockPolicy3D().getBlockCommunicator(),
          defaultMultiBlockPolicy3D().getCombinedStatistics(),
          defaultMultiBlockPolicy3D().getMultiTensorAccess<T, 3>()),
        normal((MultiBlock3D &)curvature),
        container(0)
    {
        incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif

        useSurfaceTension = !util::isZero(surfaceTension);

        freeSurfaceArgs = aggregateFreeSurfaceParams(
            lattice, rhoBar, j, mass, volumeFraction, flag, normal, helperLists, curvature,
            outsideDensity);

        initializeInterfaceLists3D<T, Descriptor>(helperLists);
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

        lattice.internalStatSubscription().subscribeSum();     // Total mass.
        lattice.internalStatSubscription().subscribeSum();     // Lost mass.
        lattice.internalStatSubscription().subscribeIntSum();  // Num interface cells.

        freeSurfaceDataProcessors();
    }

    // TODO: The default argument repelInterface is to test different
    //       methods of repelling bubbles (and droplets) from immersed boundaries.
    //       0: no repelling
    //       1: use rhoBar-j kinematic repelling
    //       2: use repelling by changing the "temporarilyProtect" flag
    //       3: use rhoBar-j kinematic repelling, but a very strong one
    //       When tests are over, this argument is meant to be disposed of.
    template <class VelFunction>
    FreeSurfaceFields3D(
        SparseBlockStructure3D const &blockStructure, Dynamics<T, Descriptor> *dynamics_,
        T rhoDefault_, T surfaceTension_, T contactAngle_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_, plint numIBIterations,
        std::vector<Array<T, 3> > const &vertices, std::vector<T> const &areas,
        std::vector<int> const &flags, VelFunction velFunction, int repelInterface = 0) :
        dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        force(force_),
        lattice(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                smallEnvelopeWidth),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiCellAccess<T, Descriptor>(), dynamics->clone()),
        helperLists(lattice),
        mass(lattice),
        flag(
            MultiBlockManagement3D(
                blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                envelopeWidthForImmersedWalls),
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
                envelopeWidthForImmersedWalls),
            defaultMultiBlockPolicy3D().getBlockCommunicator(),
            defaultMultiBlockPolicy3D().getCombinedStatistics(),
            defaultMultiBlockPolicy3D().getMultiScalarAccess<T>()),
        j(MultiBlockManagement3D(
              blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
              envelopeWidthForImmersedWalls),
          defaultMultiBlockPolicy3D().getBlockCommunicator(),
          defaultMultiBlockPolicy3D().getCombinedStatistics(),
          defaultMultiBlockPolicy3D().getMultiTensorAccess<T, 3>()),
        normal((MultiBlock3D &)curvature)
    {
        container = new MultiContainerBlock3D((MultiBlock3D &)rhoBar);
        PLB_ASSERT(container);

        incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif

        useSurfaceTension = !util::isZero(surfaceTension);

        freeSurfaceArgs = aggregateFreeSurfaceParams(
            lattice, rhoBar, j, mass, volumeFraction, flag, normal, helperLists, curvature,
            outsideDensity);

        initializeInterfaceLists3D<T, Descriptor>(helperLists);
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

        lattice.internalStatSubscription().subscribeSum();     // Total mass.
        lattice.internalStatSubscription().subscribeSum();     // Lost mass.
        lattice.internalStatSubscription().subscribeIntSum();  // Num interface cells.

        freeSurfaceDataProcessorsForImmersedWalls(
            numIBIterations, vertices, areas, flags, velFunction, repelInterface);
    }

    FreeSurfaceFields3D(FreeSurfaceFields3D<T, Descriptor> const &rhs) :
        dynamics(rhs.dynamics->clone()),
        incompressibleModel(rhs.incompressibleModel),
        rhoDefault(rhs.rhoDefault),
        surfaceTension(rhs.surfaceTension),
        contactAngle(rhs.contactAngle),
        useSurfaceTension(rhs.useSurfaceTension),
        force(rhs.force),
        lattice(rhs.lattice),
        helperLists(rhs.helperLists),
        mass(rhs.mass),
        flag(rhs.flag),
        volumeFraction(rhs.volumeFraction),
        curvature(rhs.curvature),
        outsideDensity(rhs.outsideDensity),
        rhoBar(rhs.rhoBar),
        j(rhs.j),
        normal(rhs.normal),
        rhoBarJparam(rhs.rhoBarJparam),
        freeSurfaceArgs(rhs.freeSurfaceArgs)
    {
        container = 0;
        if (rhs.container) {
            container = rhs.container->clone();
            PLB_ASSERT(container);
        }
    }

    void swap(FreeSurfaceFields3D<T, Descriptor> &rhs)
    {
        std::swap(dynamics, rhs.dynamics);
        std::swap(incompressibleModel, rhs.incompressibleModel);
        std::swap(rhoDefault, rhs.rhoDefault);
        std::swap(surfaceTension, rhs.surfaceTension);
        std::swap(contactAngle, rhs.contactAngle);
        std::swap(useSurfaceTension, rhs.useSurfaceTension);
        std::swap(force, rhs.force);
        std::swap(lattice, rhs.lattice);
        std::swap(helperLists, rhs.helperLists);
        std::swap(mass, rhs.mass);
        std::swap(flag, rhs.flag);
        std::swap(volumeFraction, rhs.volumeFraction);
        std::swap(curvature, rhs.curvature);
        std::swap(outsideDensity, rhs.outsideDensity);
        std::swap(rhoBar, rhs.rhoBar);
        std::swap(j, rhs.j);
        std::swap(normal, rhs.normal);
        std::swap(container, rhs.container);
        std::swap(rhoBarJparam, rhs.rhoBarJparam);
        std::swap(freeSurfaceArgs, rhs.freeSurfaceArgs);
    }

    FreeSurfaceFields3D<T, Descriptor> &operator=(FreeSurfaceFields3D<T, Descriptor> const &rhs)
    {
        FreeSurfaceFields3D<T, Descriptor>(rhs).swap(*this);
        return *this;
    }

    FreeSurfaceFields3D<T, Descriptor> *clone() const
    {
        return new FreeSurfaceFields3D<T, Descriptor>(*this);
    }

    ~FreeSurfaceFields3D()
    {
        delete dynamics;
        delete container;
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
        if (container) {
            container->periodicity().toggle(direction, periodic);
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
        if (container) {
            container->periodicity().toggleAll(periodic);
        }
    }

    void defaultInitialize(
        bool useConstRho = true, bool useZeroMomentum = true, bool initializeCell = true)
    {
        applyProcessingFunctional(
            new DefaultInitializeFreeSurface3D<T, Descriptor>(
                dynamics->clone(), force, rhoDefault, useConstRho, useZeroMomentum, initializeCell),
            lattice.getBoundingBox(), freeSurfaceArgs);
    }

    void partiallyDefaultInitialize(
        bool useConstRho = true, bool useZeroMomentum = true, bool initializeCell = true)
    {
        applyProcessingFunctional(
            new PartiallyDefaultInitializeFreeSurface3D<T, Descriptor>(
                dynamics->clone(), force, rhoDefault, useConstRho, useZeroMomentum, initializeCell),
            lattice.getBoundingBox(), freeSurfaceArgs);
    }

    void freeSurfaceDataProcessors()
    {
        plint pl;  // Processor level.

        /***** Initial level ******/
        pl = 0;

        integrateProcessingFunctional(
            new ExternalRhoJcollideAndStream3D<T, Descriptor>, lattice.getBoundingBox(),
            rhoBarJparam, pl);

        integrateProcessingFunctional(
            new FreeSurfaceComputeNormals3D<T, Descriptor>, lattice.getBoundingBox(),
            freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        if (useSurfaceTension) {
            integrateProcessingFunctional(
                new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngle),
                lattice.getBoundingBox(), freeSurfaceArgs, pl);

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
            new FreeSurfaceMassChange3D<T, Descriptor>, lattice.getBoundingBox(), freeSurfaceArgs,
            pl);

        integrateProcessingFunctional(
            new FreeSurfaceCompletion3D<T, Descriptor>, lattice.getBoundingBox(), freeSurfaceArgs,
            pl);

        integrateProcessingFunctional(
            new FreeSurfaceMacroscopic3D<T, Descriptor>(incompressibleModel),
            lattice.getBoundingBox(), freeSurfaceArgs, pl);

        if (useSurfaceTension) {
            integrateProcessingFunctional(
                new FreeSurfaceAddSurfaceTension3D<T, Descriptor>(
                    surfaceTension, incompressibleModel),
                lattice.getBoundingBox(), freeSurfaceArgs, pl);
        }

        integrateProcessingFunctional(
            new FreeSurfaceStabilize3D<T, Descriptor>(), lattice.getBoundingBox(), freeSurfaceArgs,
            pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), lattice.getBoundingBox(),
            freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoDefault),
            lattice.getBoundingBox(), freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(dynamics->clone(), force),
            lattice.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoDefault),
            lattice.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>(),
            lattice.getBoundingBox(), freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceComputeStatistics3D<T, Descriptor>, lattice.getBoundingBox(),
            freeSurfaceArgs, pl);

        bool useForce = !util::isZero(norm(force));
        if (useForce) {
            integrateProcessingFunctional(
                new FreeSurfaceAddExternalForce3D<T, Descriptor>(rhoDefault),
                lattice.getBoundingBox(), freeSurfaceArgs, pl);
        }
    }

    template <class VelFunction>
    void freeSurfaceDataProcessorsForImmersedWalls(
        plint numIBIterations, std::vector<Array<T, 3> > const &vertices,
        std::vector<T> const &areas, std::vector<int> const &flags, VelFunction velFunction,
        int repelInterface)
    {
        plint pl;  // Processor level.

        /***** Initial level ******/
        pl = 0;

        integrateProcessingFunctional(
            new ExternalRhoJcollideAndStream3D<T, Descriptor>, lattice.getBoundingBox(),
            rhoBarJparam, pl);

        integrateProcessingFunctional(
            new FreeSurfaceComputeNormals3D<T, Descriptor>, lattice.getBoundingBox(),
            freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        if (useSurfaceTension) {
            integrateProcessingFunctional(
                new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngle),
                lattice.getBoundingBox(), freeSurfaceArgs, pl);

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
            new FreeSurfaceMassChange3D<T, Descriptor>, lattice.getBoundingBox(), freeSurfaceArgs,
            pl);

        integrateProcessingFunctional(
            new FreeSurfaceCompletion3D<T, Descriptor>, lattice.getBoundingBox(), freeSurfaceArgs,
            pl);

        integrateProcessingFunctional(
            new FreeSurfaceMacroscopic3D<T, Descriptor>(incompressibleModel),
            lattice.getBoundingBox(), freeSurfaceArgs, pl);

        if (useSurfaceTension) {
            integrateProcessingFunctional(
                new FreeSurfaceAddSurfaceTension3D<T, Descriptor>(
                    surfaceTension, incompressibleModel),
                lattice.getBoundingBox(), freeSurfaceArgs, pl);
        }

        integrateProcessingFunctional(
            new FreeSurfaceStabilize3D<T, Descriptor>(), lattice.getBoundingBox(), freeSurfaceArgs,
            pl);

        std::vector<MultiBlock3D *> immersedWallDataArgs;
        immersedWallDataArgs.push_back(container);
        integrateProcessingFunctional(
            new InstantiateImmersedWallDataWithIndexedTagging3D<T>(vertices, areas, flags),
            container->getBoundingBox(), lattice, immersedWallDataArgs, pl);

        /***** New level ******/
        if (repelInterface == 2) {
            pl++;

            std::vector<MultiBlock3D *> tmpProtectionArgs;
            tmpProtectionArgs.push_back(&flag);
            tmpProtectionArgs.push_back(container);
            integrateProcessingFunctional(
                new TemporarilyProtectImmersedWalls3D<T>(), flag.getBoundingBox(), lattice,
                tmpProtectionArgs, pl);
        }

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), lattice.getBoundingBox(),
            freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoDefault),
            lattice.getBoundingBox(), freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(dynamics->clone(), force),
            lattice.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoDefault),
            lattice.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>(),
            lattice.getBoundingBox(), freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceComputeStatistics3D<T, Descriptor>, lattice.getBoundingBox(),
            freeSurfaceArgs, pl);

        /***** New level ******/

        for (int i = 0; i < numIBIterations; i++) {
            pl++;
            T tau = (T)1 / dynamics->getOmega();
            std::vector<MultiBlock3D *> args;
            args.push_back(&rhoBar);
            args.push_back(&j);
            args.push_back(container);
            integrateProcessingFunctional(
                new IndexedInamuroIteration3D<T, VelFunction>(
                    velFunction, tau, incompressibleModel),
                rhoBar.getBoundingBox(), lattice, args, pl);
        }

        if (repelInterface == 1 || repelInterface == 3) {
            bool strongRepelling = repelInterface == 1 ? false : true;

            // TODO: This data processor changes the momentum in the vicinity of the
            //       immersed boundary. Maybe it affects strongly the measurement of
            //       force and torque on the immersed surface. Needs to be checked.
            pl++;
            std::vector<MultiBlock3D *> args;
            args.push_back(&rhoBar);
            args.push_back(&j);
            args.push_back(&flag);
            args.push_back(container);
            integrateProcessingFunctional(
                new RepelInterfaceFromImmersedWalls3D<T, VelFunction>(
                    velFunction, rhoDefault, strongRepelling),
                rhoBar.getBoundingBox(), lattice, args, pl);
        }

        /***** New level ******/
        if (repelInterface == 2) {
            pl++;

            std::vector<MultiBlock3D *> rmProtectionArgs;
            rmProtectionArgs.push_back(&flag);
            integrateProcessingFunctional(
                new RemoveProtectionFromImmersedWalls3D<T>(), flag.getBoundingBox(), lattice,
                rmProtectionArgs, pl);
        }

        /***** New level ******/
        bool useForce = !util::isZero(norm(force));
        if (useForce) {
            pl++;

            integrateProcessingFunctional(
                new FreeSurfaceAddExternalForce3D<T, Descriptor>(rhoDefault),
                lattice.getBoundingBox(), freeSurfaceArgs, pl);
        }
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
    }

    Dynamics<T, Descriptor> *dynamics;
    bool incompressibleModel;
    T rhoDefault;
    T surfaceTension;
    T contactAngle;
    bool useSurfaceTension;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force;
    MultiBlockLattice3D<T, Descriptor> lattice;
    MultiContainerBlock3D helperLists;
    MultiScalarField3D<T> mass;
    MultiScalarField3D<int> flag;
    MultiScalarField3D<T> volumeFraction;
    MultiScalarField3D<T> curvature;
    MultiScalarField3D<T> outsideDensity;
    MultiScalarField3D<T> rhoBar;
    MultiTensorField3D<T, 3> j;
    MultiTensorField3D<T, 3> normal;
    MultiContainerBlock3D *container;
    std::vector<MultiBlock3D *> rhoBarJparam;
    std::vector<MultiBlock3D *> freeSurfaceArgs;
};

template <typename T, template <typename U> class Descriptor>
const int FreeSurfaceFields3D<T, Descriptor>::envelopeWidth =
    3;  // Necessary when we use height functions to compute the curvature,
        // or when double smoothing is used at the data processor that
        // computes the normals from the volume fraction.
// template<typename T, template<typename U> class Descriptor>
// const int FreeSurfaceFields3D<T,Descriptor>::envelopeWidth = 4; // Necessary when we use height
// functions to compute the curvature and
//  use the old contact angle algorithm.
template <typename T, template <typename U> class Descriptor>
const int FreeSurfaceFields3D<T, Descriptor>::smallEnvelopeWidth = 1;

template <typename T, template <typename U> class Descriptor>
const int FreeSurfaceFields3D<T, Descriptor>::envelopeWidthForImmersedWalls = 4;

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceSetup {
public:
    typedef T (*ContactAngleFunction)(T x, T y, T z);  // Returns the contact angle in degrees.
public:
    FreeSurfaceSetup(
        Group3D &group_, Dynamics<T, Descriptor> const &dynamics_, T rhoDefault_, T surfaceTension_,
        T contactAngle_, Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        std::string fluidname_ = "fluid1", std::string fsPrefix_ = "",
        std::string rhoBarJprefix_ = "", bool hasImmersedWalls_ = false) :
        group(group_),
        dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        contactAngleFunction(0),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_)
    {
        freeSurfaceArgs.clear();
        rhoBarJparam.clear();
        incompressibleModel = dynamics.velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif
    }

    FreeSurfaceSetup(
        Group3D &group_, Dynamics<T, Descriptor> const &dynamics_, T rhoDefault_, T surfaceTension_,
        ContactAngleFunction contactAngleFunction_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        std::string fluidname_ = "fluid1", std::string fsPrefix_ = "",
        std::string rhoBarJprefix_ = "", bool hasImmersedWalls_ = false) :
        group(group_),
        dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(-1.0),
        contactAngleFunction(contactAngleFunction_),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_)
    {
        freeSurfaceArgs.clear();
        rhoBarJparam.clear();
        incompressibleModel = dynamics.velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif
    }

    void periodicityToggle(plint direction, bool periodic)
    {
        PLB_ASSERT(direction == 0 || direction == 1 || direction == 2);

        group.get(fluidname).periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "mass").periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "flag").periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "volumeFraction").periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "curvature").periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "outsideDensity").periodicity().toggle(direction, periodic);
        group.get(rhoBarJprefix + "rhoBar").periodicity().toggle(direction, periodic);
        group.get(rhoBarJprefix + "j").periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "normal").periodicity().toggle(direction, periodic);
    }

    void periodicityToggleAll(bool periodic)
    {
        group.get(fluidname).periodicity().toggleAll(periodic);
        group.get(fsPrefix + "mass").periodicity().toggleAll(periodic);
        group.get(fsPrefix + "flag").periodicity().toggleAll(periodic);
        group.get(fsPrefix + "volumeFraction").periodicity().toggleAll(periodic);
        group.get(fsPrefix + "curvature").periodicity().toggleAll(periodic);
        group.get(fsPrefix + "outsideDensity").periodicity().toggleAll(periodic);
        group.get(rhoBarJprefix + "rhoBar").periodicity().toggleAll(periodic);
        group.get(rhoBarJprefix + "j").periodicity().toggleAll(periodic);
        group.get(fsPrefix + "normal").periodicity().toggleAll(periodic);
    }

    void createFreeSurfaceFields()
    {
        const int envelopeWidth = FreeSurfaceFields3D<T, Descriptor>::envelopeWidth;
        const int smallEnvelopeWidth = FreeSurfaceFields3D<T, Descriptor>::smallEnvelopeWidth;
        const int smallOrLargeEnvelopeWidth =
            hasImmersedWalls ? FreeSurfaceFields3D<T, Descriptor>::envelopeWidthForImmersedWalls
                             : smallEnvelopeWidth;
        const int mediumOrLargeEnvelopeWidth =
            hasImmersedWalls ? FreeSurfaceFields3D<T, Descriptor>::envelopeWidthForImmersedWalls
                             : envelopeWidth;
#ifdef PLB_DEBUG
        bool hasFluid = fieldExists(fluidname, smallEnvelopeWidth);
#endif
        PLB_ASSERT(hasFluid);
        group.get(fluidname).periodicity().toggleAll(true);
        if (!fieldExists(fsPrefix + "helperLists", 0)) {
            group.generateContainer(fsPrefix + "helperLists");
        }
        if (!fieldExists(fsPrefix + "mass", smallEnvelopeWidth)) {
            group.generateScalar<T>(fsPrefix + "mass", smallEnvelopeWidth);
        }
        group.get(fsPrefix + "mass").periodicity().toggleAll(true);
        if (!fieldExists("flag", mediumOrLargeEnvelopeWidth)) {
            group.generateScalar<int>(fsPrefix + "flag", mediumOrLargeEnvelopeWidth);
        }
        group.get(fsPrefix + "flag").periodicity().toggleAll(true);
        if (!fieldExists(fsPrefix + "volumeFraction", envelopeWidth)) {
            group.generateScalar<T>(fsPrefix + "volumeFraction", envelopeWidth);
        }
        group.get(fsPrefix + "volumeFraction").periodicity().toggleAll(true);
        if (!fieldExists(fsPrefix + "curvature", envelopeWidth)) {
            group.generateScalar<T>(fsPrefix + "curvature", envelopeWidth);
        }
        group.get(fsPrefix + "curvature").periodicity().toggleAll(true);
        if (!fieldExists(fsPrefix + "outsideDensity", envelopeWidth)) {
            group.generateScalar<T>(fsPrefix + "outsideDensity", envelopeWidth);
        }
        group.get(fsPrefix + "outsideDensity").periodicity().toggleAll(true);
        if (!fieldExists(rhoBarJprefix + "rhoBar", smallOrLargeEnvelopeWidth)) {
            group.generateScalar<T>(rhoBarJprefix + "rhoBar", smallOrLargeEnvelopeWidth);
        }
        group.get(rhoBarJprefix + "rhoBar").periodicity().toggleAll(true);
        if (!fieldExists(rhoBarJprefix + "j", smallOrLargeEnvelopeWidth)) {
            group.generateTensor<T, 3>(rhoBarJprefix + "j", smallOrLargeEnvelopeWidth);
        }
        group.get(rhoBarJprefix + "j").periodicity().toggleAll(true);
        if (!fieldExists(fsPrefix + "normal", envelopeWidth)) {
            group.generateTensor<T, 3>(fsPrefix + "normal", envelopeWidth);
        }
        group.get(fsPrefix + "normal").periodicity().toggleAll(true);
        if (hasImmersedWalls) {
            if (!fieldExists(fsPrefix + "ibm_container", smallOrLargeEnvelopeWidth)) {
                group.generateContainer(fsPrefix + "ibm_container", smallOrLargeEnvelopeWidth);
            }
        }

        freeSurfaceArgs.push_back(&group.get(fluidname));
        freeSurfaceArgs.push_back(&group.get(rhoBarJprefix + "rhoBar"));
        freeSurfaceArgs.push_back(&group.get(rhoBarJprefix + "j"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "mass"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "volumeFraction"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "flag"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "normal"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "helperLists"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "curvature"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "outsideDensity"));

        rhoBarJparam.push_back(&group.get(fluidname));
        rhoBarJparam.push_back(&group.get(rhoBarJprefix + "rhoBar"));
        rhoBarJparam.push_back(&group.get(rhoBarJprefix + "j"));
    }

    Actions3D freeSurfaceActions(
        std::map<std::string, plint> *ids = 0, T characteristicLength = (T)0, T criticalWe = (T)0,
        T criticalKn = (T)0)
    {
        Actions3D actions;

        bool useSurfaceTension = !util::isZero(surfaceTension);
        bool useWeberModel = !util::isZero(characteristicLength) && !util::isZero(criticalWe)
                             && !util::isZero(criticalKn);

        plint id = -1;
        id = actions.addBlock(group.get(fluidname));  // 0
        if (ids != 0)
            (*ids)[fluidname] = id;
        id = actions.addBlock(group.get(rhoBarJprefix + "rhoBar"));  // 1
        if (ids != 0)
            (*ids)[rhoBarJprefix + "rhoBar"] = id;
        id = actions.addBlock(group.get(rhoBarJprefix + "j"));  // 2
        if (ids != 0)
            (*ids)[rhoBarJprefix + "j"] = id;
        id = actions.addBlock(group.get(fsPrefix + "mass"));  // 3
        if (ids != 0)
            (*ids)[fsPrefix + "mass"] = id;
        id = actions.addBlock(group.get(fsPrefix + "volumeFraction"));  // 4
        if (ids != 0)
            (*ids)[fsPrefix + "volumeFraction"] = id;
        id = actions.addBlock(group.get(fsPrefix + "flag"));  // 5
        if (ids != 0)
            (*ids)[fsPrefix + "flag"] = id;
        id = actions.addBlock(group.get(fsPrefix + "normal"));  // 6
        if (ids != 0)
            (*ids)[fsPrefix + "normal"] = id;
        id = actions.addBlock(group.get(fsPrefix + "helperLists"));  // 7
        if (ids != 0)
            (*ids)[fsPrefix + "helperLists"] = id;
        id = actions.addBlock(group.get(fsPrefix + "curvature"));  // 8
        if (ids != 0)
            (*ids)[fsPrefix + "curvature"] = id;
        id = actions.addBlock(group.get(fsPrefix + "outsideDensity"));  // 9
        if (ids != 0)
            (*ids)[fsPrefix + "outsideDensity"] = id;

        std::vector<plint> freeSurfaceBlocks;
        for (plint i = 0; i < 10; ++i)
            freeSurfaceBlocks.push_back(i);

        initializeInterfaceLists3D<T, Descriptor>(group.getContainer(fsPrefix + "helperLists"));
        // setToConstant(group.getScalar<int>("flag"), group.getBoundingBox(),
        // (int)freeSurfaceFlag::empty); setToConstant(group.getScalar<T>("outsideDensity"),
        // group.getBoundingBox(), rhoDefault);

        std::vector<plint> rhoBarJblocks;
        for (plint i = 0; i < 3; ++i)
            rhoBarJblocks.push_back(i);

        group.get(fluidname).internalStatSubscription().subscribeSum();     // Total mass.
        group.get(fluidname).internalStatSubscription().subscribeSum();     // Lost mass.
        group.get(fluidname).internalStatSubscription().subscribeIntSum();  // Num interface cells.

        actions.addProcessor(
            new ExternalRhoJcollideAndStream3D<T, Descriptor>, rhoBarJblocks,
            group.getBoundingBox());
        actions.addProcessor(
            new FreeSurfaceComputeNormals3D<T, Descriptor>, freeSurfaceBlocks,
            group.getBoundingBox());
        actions.addCommunication(0, modif::staticVariables);  // fluid1
        actions.addCommunication(6, modif::staticVariables);  // normal

        if (useSurfaceTension) {
            if (contactAngleFunction == 0) {
                actions.addProcessor(
                    new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngle),
                    freeSurfaceBlocks, group.getBoundingBox());
            } else {
                actions.addProcessor(
                    new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngleFunction),
                    freeSurfaceBlocks, group.getBoundingBox());
            }
        }
        actions.addProcessor(
            new FreeSurfaceMassChange3D<T, Descriptor>, freeSurfaceBlocks, group.getBoundingBox());
        actions.addProcessor(
            new FreeSurfaceCompletion3D<T, Descriptor>, freeSurfaceBlocks, group.getBoundingBox());
        actions.addProcessor(
            new FreeSurfaceMacroscopic3D<T, Descriptor>(incompressibleModel), freeSurfaceBlocks,
            group.getBoundingBox());
        if (useSurfaceTension) {
            if (!useWeberModel) {
                actions.addProcessor(
                    new FreeSurfaceAddSurfaceTension3D<T, Descriptor>(
                        surfaceTension, incompressibleModel),
                    freeSurfaceBlocks, group.getBoundingBox());
            } else {
                actions.addProcessor(
                    new FreeSurfaceAddSurfaceTensionWeberModel3D<T, Descriptor>(
                        surfaceTension, rhoDefault, characteristicLength, criticalWe, criticalKn,
                        incompressibleModel),
                    freeSurfaceBlocks, group.getBoundingBox());
            }
        }
        actions.addProcessor(
            new FreeSurfaceStabilize3D<T, Descriptor>(), freeSurfaceBlocks, group.getBoundingBox());
        actions.addCommunication(0, modif::staticVariables);  // fluid1
        actions.addCommunication(1, modif::staticVariables);  // rhoBar.
        actions.addCommunication(2, modif::staticVariables);  // j.
        actions.addCommunication(3, modif::staticVariables);  // Mass.
        actions.addCommunication(4, modif::staticVariables);  // Volume fraction.
        if (useSurfaceTension) {
            actions.addCommunication(8, modif::staticVariables);  // Curvature.
        }

        actions.addProcessor(
            new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), freeSurfaceBlocks,
            group.getBoundingBox());
        actions.addProcessor(
            new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoDefault), freeSurfaceBlocks,
            group.getBoundingBox());
        actions.addProcessor(
            new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(dynamics.clone(), force),
            freeSurfaceBlocks, group.getBoundingBox());

        actions.addCommunication(1, modif::staticVariables);  // rhoBar.
        actions.addCommunication(3, modif::staticVariables);  // Mass.
        actions.addCommunication(5, modif::staticVariables);  // Flag-status.
        actions.addCommunication(7, modif::staticVariables);  // Interface-lists.

        actions.addProcessor(
            new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoDefault),
            freeSurfaceBlocks, group.getBoundingBox());
        actions.addCommunication(1, modif::staticVariables);  // rhoBar.
        actions.addCommunication(3, modif::staticVariables);  // Mass.
        actions.addCommunication(5, modif::staticVariables);  // Flag-status.

        actions.addProcessor(
            new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>, freeSurfaceBlocks,
            group.getBoundingBox());
        actions.addProcessor(
            new FreeSurfaceComputeStatistics3D<T, Descriptor>, freeSurfaceBlocks,
            group.getBoundingBox());
        bool useForce = !util::isZero(norm(force));
        if (useForce) {
            actions.addProcessor(
                new FreeSurfaceAddExternalForce3D<T, Descriptor>(rhoDefault), freeSurfaceBlocks,
                group.getBoundingBox());
        }
        actions.addCommunication(0, modif::dataStructure);    // Fluid.
        actions.addCommunication(1, modif::staticVariables);  // rhoBar.
        actions.addCommunication(2, modif::staticVariables);  // j.
        actions.addCommunication(3, modif::staticVariables);  // Mass.
        actions.addCommunication(4, modif::staticVariables);  // Volume fraction.

        actions.addEvaluateStats(0);                 // Fluid.
        actions.addIncrementTime<T, Descriptor>(0);  // Fluid.

        return actions;
    }

    template <class VelFunction>
    Actions3D immersedWallActions(  // TODO: maybe remove container?
        plint numIBIterations, [[maybe_unused]] MultiContainerBlock3D &container,
        std::vector<Array<T, 3> > const &vertices, std::vector<T> const &areas,
        std::vector<int> const &flags, VelFunction velFunction, bool strongRepelling)
    {
        Actions3D actions;

        plint containerID = actions.addBlock(group.get(fsPrefix + "ibm_container"));
        plint rhoBarID = actions.addBlock(group.get(rhoBarJprefix + "rhoBar"));
        plint jID = actions.addBlock(group.get(rhoBarJprefix + "j"));
        plint flagID = actions.addBlock(group.get(fsPrefix + "flag"));

        actions.addProcessor(
            new InstantiateImmersedWallDataWithIndexedTagging3D<T>(vertices, areas, flags),
            containerID, group.getBoundingBox());
        T tau = (T)1 / dynamics.getOmega();
        for (int i = 0; i < numIBIterations; i++) {
            actions.addProcessor(
                new IndexedInamuroIteration3D<T, VelFunction>(
                    velFunction, tau, incompressibleModel),
                rhoBarID, jID, containerID, group.getBoundingBox());
            actions.addCommunication(jID, modif::staticVariables);
        }
        actions.addProcessor(
            new RepelInterfaceFromImmersedWalls3D<T, VelFunction>(
                velFunction, rhoDefault, strongRepelling),
            rhoBarID, jID, flagID, containerID, group.getBoundingBox());
        actions.addCommunication(rhoBarID, modif::staticVariables);
        actions.addCommunication(jID, modif::staticVariables);

        return actions;
    }

    // If surfaceTensionField != 0, then the FreeSurfaceAddSurfaceTensionFromScalarField3D processor
    // is integrated instead of the FreeSurfaceAddSurfaceTension3D one. Contrary to the Palabos
    // convention, this function does NOT take ownership of the surfaceTensionField despite the fact
    // that a pointer is passed. The caller is responsible for its memory management.
    void createFreeSurfaceProcessors(
        plint initialProcessorLevel = 0, MultiScalarField3D<T> *surfaceTensionField = 0)
    {
        PLB_ASSERT(initialProcessorLevel >= 0);

        bool useSurfaceTension = (surfaceTensionField != 0 || !util::isZero(surfaceTension));

        initializeInterfaceLists3D<T, Descriptor>(group.getContainer(fsPrefix + "helperLists"));
        // setToConstant(group.getScalar<int>("flag"), group.getBoundingBox(),
        // (int)freeSurfaceFlag::empty); setToConstant(group.getScalar<T>("outsideDensity"),
        // group.getBoundingBox(), rhoDefault);

        group.get(fluidname).internalStatSubscription().subscribeSum();     // Total mass.
        group.get(fluidname).internalStatSubscription().subscribeSum();     // Lost mass.
        group.get(fluidname).internalStatSubscription().subscribeIntSum();  // Num interface cells.

        plint pl = initialProcessorLevel;  // Processor level.

        /***** Initial level ******/
        if (pl == 0) {
            integrateProcessingFunctional(
                new ExternalRhoJcollideAndStream3D<T, Descriptor>, group.getBoundingBox(),
                rhoBarJparam, pl);
        }

        integrateProcessingFunctional(
            new FreeSurfaceComputeNormals3D<T, Descriptor>, group.getBoundingBox(), freeSurfaceArgs,
            pl);

        /***** New level ******/
        pl++;

        if (useSurfaceTension) {
            if (contactAngleFunction == 0) {
                integrateProcessingFunctional(
                    new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngle),
                    group.getBoundingBox(), freeSurfaceArgs, pl);
            } else {
                integrateProcessingFunctional(
                    new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngleFunction),
                    group.getBoundingBox(), freeSurfaceArgs, pl);
            }

            // To change to the curvature calculation with height functions, uncomment the next data
            // processor and comment out the two previous ones. If only the next data processor is
            // used and there is no surface tension, the normals are not computed at all. Be careful
            // if you intent to use the normals and do not have the surface tension algorithm
            // enabled.
            // integrateProcessingFunctional (
            //        new FreeSurfaceGeometry3D<T,Descriptor>(contactAngle),
            //        group.getBoundingBox(), freeSurfaceArgs, pl );
        }

        integrateProcessingFunctional(
            new FreeSurfaceMassChange3D<T, Descriptor>, group.getBoundingBox(), freeSurfaceArgs,
            pl);

        integrateProcessingFunctional(
            new FreeSurfaceCompletion3D<T, Descriptor>, group.getBoundingBox(), freeSurfaceArgs,
            pl);

        integrateProcessingFunctional(
            new FreeSurfaceMacroscopic3D<T, Descriptor>(incompressibleModel),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        if (useSurfaceTension) {
            if (surfaceTensionField == 0) {
                integrateProcessingFunctional(
                    new FreeSurfaceAddSurfaceTension3D<T, Descriptor>(
                        surfaceTension, incompressibleModel),
                    group.getBoundingBox(), freeSurfaceArgs, pl);
            } else {
                std::vector<MultiBlock3D *> args(freeSurfaceArgs);
                args.push_back(surfaceTensionField);

                integrateProcessingFunctional(
                    new FreeSurfaceAddSurfaceTensionFromScalarField3D<T, Descriptor>(
                        incompressibleModel),
                    group.getBoundingBox(), args, pl);
            }
        }

        integrateProcessingFunctional(
            new FreeSurfaceStabilize3D<T, Descriptor>(), group.getBoundingBox(), freeSurfaceArgs,
            pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), group.getBoundingBox(),
            freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(rhoDefault),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(dynamics.clone(), force),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(rhoDefault),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>(), group.getBoundingBox(),
            freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceComputeStatistics3D<T, Descriptor>, group.getBoundingBox(),
            freeSurfaceArgs, pl);

        bool useForce = !util::isZero(norm(force));
        if (useForce) {
            integrateProcessingFunctional(
                new FreeSurfaceAddExternalForce3D<T, Descriptor>(rhoDefault),
                group.getBoundingBox(), freeSurfaceArgs, pl);
        }
    }

    void defaultInitialize(
        bool useConstRho = true, bool useZeroMomentum = true, bool initializeCell = true)
    {
        applyProcessingFunctional(
            new DefaultInitializeFreeSurface3D<T, Descriptor>(
                dynamics.clone(), force, rhoDefault, useConstRho, useZeroMomentum, initializeCell),
            group.get(fluidname).getBoundingBox(), freeSurfaceArgs);
    }

    void partiallyDefaultInitialize(
        bool useConstRho = true, bool useZeroMomentum = true, bool initializeCell = true)
    {
        applyProcessingFunctional(
            new PartiallyDefaultInitializeFreeSurface3D<T, Descriptor>(
                dynamics.clone(), force, rhoDefault, useConstRho, useZeroMomentum, initializeCell),
            group.get(fluidname).getBoundingBox(), freeSurfaceArgs);
    }

    Group3D &getGroup()
    {
        return group;
    }

    Group3D const &getGroup() const
    {
        return group;
    }

    Dynamics<T, Descriptor> const &getDynamics() const
    {
        return dynamics;
    }

    T getRhoDefault() const
    {
        return rhoDefault;
    }

    T getSurfaceTension() const
    {
        return surfaceTension;
    }

    T getContactAngle() const
    {
        return contactAngle;
    }

    ContactAngleFunction getContactAngleFunction() const
    {
        return contactAngleFunction;
    }

    Array<T, Descriptor<T>::ExternalField::sizeOfForce> const &getForce() const
    {
        return force;
    }

    bool getIncompressibleModel() const
    {
        return incompressibleModel;
    }

private:
    // TODO: maybe replace ASSERTS with exceptions.
    bool fieldExists(std::string name, [[maybe_unused]] plint envelopeWidth)
    {
        if (group.hasBlock(name)) {
            PLB_ASSERT(
                group.get(name).getMultiBlockManagement().getEnvelopeWidth() >= envelopeWidth);
            return true;
        }
        return false;
    }

private:
    Group3D &group;
    Dynamics<T, Descriptor> const &dynamics;
    bool incompressibleModel;
    T rhoDefault, surfaceTension, contactAngle;
    ContactAngleFunction contactAngleFunction;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force;
    std::string fluidname, fsPrefix, rhoBarJprefix;
    bool hasImmersedWalls;

public:
    std::vector<MultiBlock3D *> freeSurfaceArgs;
    std::vector<MultiBlock3D *> rhoBarJparam;
};

template <typename T, template <typename U> class Descriptor>
class FreeSurfaceWrapper {
public:
    typedef T (*ContactAngleFunction)(T x, T y, T z);  // Returns the contact angle in degrees.

public:
    FreeSurfaceWrapper(
        Group3D &group_, Dynamics<T, Descriptor> *dynamics_,
        Dynamics<T, Descriptor> *emptyNodeDynamics_, T rhoDefault_, T surfaceTension_,
        T contactAngle_, Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        std::string fluidname_ = "fluid1", std::string fsPrefix_ = "",
        std::string rhoBarJprefix_ = "", bool hasImmersedWalls_ = false) :
        group(group_),
        lattice(0),
        dynamics(dynamics_),
        emptyNodeDynamics(emptyNodeDynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        contactAngleFunction(0),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_)
    {
        freeSurfaceArgs.clear();
        rhoBarJparam.clear();
        incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif
    }

    FreeSurfaceWrapper(
        Group3D &group_, Dynamics<T, Descriptor> *dynamics_, T rhoDefault_, T surfaceTension_,
        T contactAngle_, Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        std::string fluidname_ = "fluid1", std::string fsPrefix_ = "",
        std::string rhoBarJprefix_ = "", bool hasImmersedWalls_ = false) :
        group(group_),
        lattice(0),
        dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        contactAngleFunction(0),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_)
    {
        emptyNodeDynamics = new NoDynamics<T, Descriptor>(rhoDefault);
        freeSurfaceArgs.clear();
        rhoBarJparam.clear();
        incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif
    }

    FreeSurfaceWrapper(
        Group3D &group_, Dynamics<T, Descriptor> *dynamics_,
        Dynamics<T, Descriptor> *emptyNodeDynamics_, T rhoDefault_, T surfaceTension_,
        ContactAngleFunction contactAngleFunction_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        std::string fluidname_ = "fluid1", std::string fsPrefix_ = "",
        std::string rhoBarJprefix_ = "", bool hasImmersedWalls_ = false) :
        group(group_),
        lattice(0),
        dynamics(dynamics_),
        emptyNodeDynamics(emptyNodeDynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(-1.0),
        contactAngleFunction(contactAngleFunction_),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_)
    {
        freeSurfaceArgs.clear();
        rhoBarJparam.clear();
        incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif
    }

    FreeSurfaceWrapper(
        Group3D &group_, Dynamics<T, Descriptor> *dynamics_, T rhoDefault_, T surfaceTension_,
        ContactAngleFunction contactAngleFunction_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        std::string fluidname_ = "fluid1", std::string fsPrefix_ = "",
        std::string rhoBarJprefix_ = "", bool hasImmersedWalls_ = false) :
        group(group_),
        lattice(0),
        dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(-1.0),
        contactAngleFunction(contactAngleFunction_),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_)
    {
        emptyNodeDynamics = new NoDynamics<T, Descriptor>(rhoDefault);
        freeSurfaceArgs.clear();
        rhoBarJparam.clear();
        incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif
    }

    // With this constructor, the fluid lattice is provided separately and is not contained in the
    // group.
    FreeSurfaceWrapper(
        Group3D &group_, MultiBlockLattice3D<T, Descriptor> &lattice_,
        Dynamics<T, Descriptor> *dynamics_, Dynamics<T, Descriptor> *emptyNodeDynamics_,
        T rhoDefault_, T surfaceTension_, T contactAngle_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        std::string fluidname_ = "fluid1", std::string fsPrefix_ = "",
        std::string rhoBarJprefix_ = "", bool hasImmersedWalls_ = false) :
        group(group_),
        lattice(&lattice_),
        dynamics(dynamics_),
        emptyNodeDynamics(emptyNodeDynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        contactAngleFunction(0),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_)
    {
        freeSurfaceArgs.clear();
        rhoBarJparam.clear();
        incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif
    }

    // With this constructor, the fluid lattice is provided separately and is not contained in the
    // group.
    FreeSurfaceWrapper(
        Group3D &group_, MultiBlockLattice3D<T, Descriptor> &lattice_,
        Dynamics<T, Descriptor> *dynamics_, T rhoDefault_, T surfaceTension_, T contactAngle_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        std::string fluidname_ = "fluid1", std::string fsPrefix_ = "",
        std::string rhoBarJprefix_ = "", bool hasImmersedWalls_ = false) :
        group(group_),
        lattice(&lattice_),
        dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(contactAngle_),
        contactAngleFunction(0),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_)
    {
        emptyNodeDynamics = new NoDynamics<T, Descriptor>(rhoDefault);
        freeSurfaceArgs.clear();
        rhoBarJparam.clear();
        incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif
    }

    // With this constructor, the fluid lattice is provided separately and is not contained in the
    // group.
    FreeSurfaceWrapper(
        Group3D &group_, MultiBlockLattice3D<T, Descriptor> &lattice_,
        Dynamics<T, Descriptor> *dynamics_, Dynamics<T, Descriptor> *emptyNodeDynamics_,
        T rhoDefault_, T surfaceTension_, ContactAngleFunction contactAngleFunction_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        std::string fluidname_ = "fluid1", std::string fsPrefix_ = "",
        std::string rhoBarJprefix_ = "", bool hasImmersedWalls_ = false) :
        group(group_),
        lattice(&lattice_),
        dynamics(dynamics_),
        emptyNodeDynamics(emptyNodeDynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(-1.0),
        contactAngleFunction(contactAngleFunction_),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_)
    {
        freeSurfaceArgs.clear();
        rhoBarJparam.clear();
        incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif
    }

    // With this constructor, the fluid lattice is provided separately and is not contained in the
    // group.
    FreeSurfaceWrapper(
        Group3D &group_, MultiBlockLattice3D<T, Descriptor> &lattice_,
        Dynamics<T, Descriptor> *dynamics_, T rhoDefault_, T surfaceTension_,
        ContactAngleFunction contactAngleFunction_,
        Array<T, Descriptor<T>::ExternalField::sizeOfForce> force_,
        std::string fluidname_ = "fluid1", std::string fsPrefix_ = "",
        std::string rhoBarJprefix_ = "", bool hasImmersedWalls_ = false) :
        group(group_),
        lattice(&lattice_),
        dynamics(dynamics_),
        rhoDefault(rhoDefault_),
        surfaceTension(surfaceTension_),
        contactAngle(-1.0),
        contactAngleFunction(contactAngleFunction_),
        force(force_),
        fluidname(fluidname_),
        fsPrefix(fsPrefix_),
        rhoBarJprefix(rhoBarJprefix_),
        hasImmersedWalls(hasImmersedWalls_)
    {
        emptyNodeDynamics = new NoDynamics<T, Descriptor>(rhoDefault);
        freeSurfaceArgs.clear();
        rhoBarJparam.clear();
        incompressibleModel = dynamics->velIsJ();
#ifdef PLB_DEBUG
        if (incompressibleModel) {
            // Incompressible: rho0=1
            PLB_ASSERT(util::isOne(rhoDefault));
        }
#endif
    }

    FreeSurfaceWrapper(FreeSurfaceWrapper<T, Descriptor> const &rhs) :
        group(rhs.group),
        lattice(rhs.lattice),
        dynamics(rhs.dynamics->clone()),
        emptyNodeDynamics(rhs.emptyNodeDynamics->clone()),
        incompressibleModel(rhs.incompressibleModel),
        rhoDefault(rhs.rhoDefault),
        surfaceTension(rhs.surfaceTension),
        contactAngle(rhs.contactAngle),
        contactAngleFunction(rhs.contactAngleFunction),
        force(rhs.force),
        fluidname(rhs.fluidname),
        fsPrefix(rhs.fsPrefix),
        rhoBarJprefix(rhs.rhoBarJprefix),
        hasImmersedWalls(rhs.hasImmersedWalls),
        reductionData(rhs.reductionData),
#ifdef PLB_MPI_PARALLEL
        reductionCommunicator(rhs.reductionCommunicator),
#endif
        freeSurfaceArgs(rhs.freeSurfaceArgs),
        rhoBarJparam(rhs.rhoBarJparam)
    { }

    void swap(FreeSurfaceWrapper<T, Descriptor> &rhs)
    {
        std::swap(group, rhs.group);
        std::swap(lattice, rhs.lattice);
        std::swap(dynamics, rhs.dynamics);
        std::swap(emptyNodeDynamics, rhs.emptyNodeDynamics);
        std::swap(incompressibleModel, rhs.incompressibleModel);
        std::swap(rhoDefault, rhs.rhoDefault);
        std::swap(surfaceTension, rhs.surfaceTension);
        std::swap(contactAngle, rhs.contactAngle);
        std::swap(contactAngleFunction, rhs.contactAngleFunction);
        std::swap(force, rhs.force);
        std::swap(fluidname, rhs.fluidname);
        std::swap(fsPrefix, rhs.fsPrefix);
        std::swap(rhoBarJprefix, rhs.rhoBarJprefix);
        std::swap(hasImmersedWalls, rhs.hasImmersedWalls);
        std::swap(reductionData, rhs.reductionData);
#ifdef PLB_MPI_PARALLEL
        std::swap(reductionCommunicator, rhs.reductionCommunicator);
#endif
        std::swap(freeSurfaceArgs, rhs.freeSurfaceArgs);
        std::swap(rhoBarJparam, rhs.rhoBarJparam);
    }

    FreeSurfaceWrapper<T, Descriptor> &operator=(FreeSurfaceWrapper<T, Descriptor> const &rhs)
    {
        FreeSurfaceWrapper<T, Descriptor>(rhs).swap(*this);
        return *this;
    }

    FreeSurfaceWrapper<T, Descriptor> *clone() const
    {
        return new FreeSurfaceWrapper<T, Descriptor>(*this);
    }

    ~FreeSurfaceWrapper()
    {
        delete dynamics;
        delete emptyNodeDynamics;
    }

    void periodicityToggle(plint direction, bool periodic)
    {
        PLB_ASSERT(direction == 0 || direction == 1 || direction == 2);

        (lattice ? *lattice : group.get(fluidname)).periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "mass").periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "flag").periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "volumeFraction").periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "curvature").periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "outsideDensity").periodicity().toggle(direction, periodic);
        group.get(rhoBarJprefix + "rhoBar").periodicity().toggle(direction, periodic);
        group.get(rhoBarJprefix + "j").periodicity().toggle(direction, periodic);
        group.get(fsPrefix + "normal").periodicity().toggle(direction, periodic);
    }

    void periodicityToggleAll(bool periodic)
    {
        (lattice ? *lattice : group.get(fluidname)).periodicity().toggleAll(periodic);
        group.get(fsPrefix + "mass").periodicity().toggleAll(periodic);
        group.get(fsPrefix + "flag").periodicity().toggleAll(periodic);
        group.get(fsPrefix + "volumeFraction").periodicity().toggleAll(periodic);
        group.get(fsPrefix + "curvature").periodicity().toggleAll(periodic);
        group.get(fsPrefix + "outsideDensity").periodicity().toggleAll(periodic);
        group.get(rhoBarJprefix + "rhoBar").periodicity().toggleAll(periodic);
        group.get(rhoBarJprefix + "j").periodicity().toggleAll(periodic);
        group.get(fsPrefix + "normal").periodicity().toggleAll(periodic);
    }

    void createFreeSurfaceFields()
    {
        const int envelopeWidth = FreeSurfaceFields3D<T, Descriptor>::envelopeWidth;
        const int smallEnvelopeWidth = FreeSurfaceFields3D<T, Descriptor>::smallEnvelopeWidth;
        const int smallOrLargeEnvelopeWidth =
            hasImmersedWalls ? FreeSurfaceFields3D<T, Descriptor>::envelopeWidthForImmersedWalls
                             : smallEnvelopeWidth;
        const int mediumOrLargeEnvelopeWidth =
            hasImmersedWalls ? FreeSurfaceFields3D<T, Descriptor>::envelopeWidthForImmersedWalls
                             : envelopeWidth;
#ifdef PLB_DEBUG
        bool hasFluid = fieldExists(fluidname, smallEnvelopeWidth);
#endif
        PLB_ASSERT((hasFluid && !lattice) || (!hasFluid && lattice));
        (lattice ? *lattice : group.get(fluidname)).periodicity().toggleAll(true);

        plint gridLevel =
            (lattice ? lattice->getMultiBlockManagement().getRefinementLevel()
                     : group.getMultiBlockManagement().getRefinementLevel());

        if (group.getNumBlocks()) {
            if (!fieldExists(fsPrefix + "mass", smallEnvelopeWidth)) {
                group.generateScalar<T>(fsPrefix + "mass", smallEnvelopeWidth, gridLevel);
            }
            group.get(fsPrefix + "mass").periodicity().toggleAll(true);
        } else {
            MultiScalarField3D<T> *mass =
                generateMultiScalarField<T>((MultiBlock3D &)*lattice, smallEnvelopeWidth).release();
            group.add(mass, fsPrefix + "mass");
            group.get(fsPrefix + "mass").periodicity().toggleAll(true);
        }

        if (!fieldExists(fsPrefix + "helperLists", 0)) {
            group.generateContainer(fsPrefix + "helperLists", smallEnvelopeWidth, gridLevel);
        }
        if (!fieldExists("flag", mediumOrLargeEnvelopeWidth)) {
            group.generateScalar<int>(fsPrefix + "flag", mediumOrLargeEnvelopeWidth, gridLevel);
        }
        group.get(fsPrefix + "flag").periodicity().toggleAll(true);
        if (!fieldExists(fsPrefix + "volumeFraction", envelopeWidth)) {
            group.generateScalar<T>(fsPrefix + "volumeFraction", envelopeWidth, gridLevel);
        }
        group.get(fsPrefix + "volumeFraction").periodicity().toggleAll(true);
        if (!fieldExists(fsPrefix + "curvature", envelopeWidth)) {
            group.generateScalar<T>(fsPrefix + "curvature", envelopeWidth, gridLevel);
        }
        group.get(fsPrefix + "curvature").periodicity().toggleAll(true);
        if (!fieldExists(fsPrefix + "outsideDensity", envelopeWidth)) {
            group.generateScalar<T>(fsPrefix + "outsideDensity", envelopeWidth, gridLevel);
        }
        group.get(fsPrefix + "outsideDensity").periodicity().toggleAll(true);
        if (!fieldExists(rhoBarJprefix + "rhoBar", smallOrLargeEnvelopeWidth)) {
            group.generateScalar<T>(rhoBarJprefix + "rhoBar", smallOrLargeEnvelopeWidth, gridLevel);
        }
        group.get(rhoBarJprefix + "rhoBar").periodicity().toggleAll(true);
        if (!fieldExists(rhoBarJprefix + "j", smallOrLargeEnvelopeWidth)) {
            group.generateTensor<T, 3>(rhoBarJprefix + "j", smallOrLargeEnvelopeWidth, gridLevel);
        }
        group.get(rhoBarJprefix + "j").periodicity().toggleAll(true);
        if (!fieldExists(fsPrefix + "normal", envelopeWidth)) {
            group.generateTensor<T, 3>(fsPrefix + "normal", envelopeWidth, gridLevel);
        }
        group.get(fsPrefix + "normal").periodicity().toggleAll(true);
        if (hasImmersedWalls) {
            if (!fieldExists(fsPrefix + "ibm_container", smallOrLargeEnvelopeWidth)) {
                group.generateContainer(
                    fsPrefix + "ibm_container", smallOrLargeEnvelopeWidth, gridLevel);
            }
        }

        freeSurfaceArgs.push_back(lattice ? lattice : &group.get(fluidname));
        freeSurfaceArgs.push_back(&group.get(rhoBarJprefix + "rhoBar"));
        freeSurfaceArgs.push_back(&group.get(rhoBarJprefix + "j"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "mass"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "volumeFraction"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "flag"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "normal"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "helperLists"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "curvature"));
        freeSurfaceArgs.push_back(&group.get(fsPrefix + "outsideDensity"));

        rhoBarJparam.push_back(lattice ? lattice : &group.get(fluidname));
        rhoBarJparam.push_back(&group.get(rhoBarJprefix + "rhoBar"));
        rhoBarJparam.push_back(&group.get(rhoBarJprefix + "j"));
    }

    // If surfaceTensionField != 0, then the FreeSurfaceAddSurfaceTensionFromScalarField3D processor
    // is integrated instead of the FreeSurfaceAddSurfaceTension3D one. Contrary to the Palabos
    // convention, this function does NOT take ownership of the surfaceTensionField despite the fact
    // that a pointer is passed. The caller is responsible for its memory management.
    void createFreeSurfaceProcessors(
        plint initialProcessorLevel, bool collideAndStreamAtTheEndOfCycle,
        BoxProcessingFunctional3D *forceProcessor, std::vector<MultiBlock3D *> blocks,
        MultiScalarField3D<T> *surfaceTensionField)
    {
        PLB_ASSERT(initialProcessorLevel >= 0);

        bool useForce = !util::isZero(norm(force));
        PLB_ASSERT(!(useForce && forceProcessor));

        bool useSurfaceTension = (surfaceTensionField != 0 || !util::isZero(surfaceTension));

        initializeInterfaceLists3D<T, Descriptor>(group.getContainer(fsPrefix + "helperLists"));
        // setToConstant(group.getScalar<int>("flag"), group.getBoundingBox(),
        // (int)freeSurfaceFlag::empty); setToConstant(group.getScalar<T>("outsideDensity"),
        // group.getBoundingBox(), rhoDefault);

        MultiBlockLattice3D<T, Descriptor> &fluidLattice =
            (lattice ? *lattice : group.getLattice<T, Descriptor>(fluidname));

        plint pl = initialProcessorLevel;  // Processor level.

        /***** Initial level ******/
        if (pl == 0 && !collideAndStreamAtTheEndOfCycle) {
            integrateProcessingFunctional(
                new ExternalRhoJcollideAndStream3D<T, Descriptor>, group.getBoundingBox(),
                rhoBarJparam, pl);
        }

        integrateProcessingFunctional(
            new FreeSurfaceComputeNormals3D<T, Descriptor>, group.getBoundingBox(), freeSurfaceArgs,
            pl);

        /***** New level ******/
        pl++;

        if (useSurfaceTension) {
            if (contactAngleFunction == 0) {
                integrateProcessingFunctional(
                    new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngle),
                    group.getBoundingBox(), freeSurfaceArgs, pl);
            } else {
                integrateProcessingFunctional(
                    new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngleFunction),
                    group.getBoundingBox(), freeSurfaceArgs, pl);
            }

            // To change to the curvature calculation with height functions, uncomment the next data
            // processor and comment out the two previous ones. If only the next data processor is
            // used and there is no surface tension, the normals are not computed at all. Be careful
            // if you intent to use the normals and do not have the surface tension algorithm
            // enabled.
            // integrateProcessingFunctional (
            //        new FreeSurfaceGeometry3D<T,Descriptor>(contactAngle),
            //        group.getBoundingBox(), freeSurfaceArgs, pl );
        }

        integrateProcessingFunctional(
            new FreeSurfaceMassChange3D<T, Descriptor>, group.getBoundingBox(), freeSurfaceArgs,
            pl);

        integrateProcessingFunctional(
            new FreeSurfaceCompletion3D<T, Descriptor>, group.getBoundingBox(), freeSurfaceArgs,
            pl);

        integrateProcessingFunctional(
            new FreeSurfaceMacroscopicWithoutLostMassReDistribution3D<T, Descriptor>(
                incompressibleModel),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        if (useSurfaceTension) {
            if (surfaceTensionField == 0) {
                integrateProcessingFunctional(
                    new FreeSurfaceAddSurfaceTension3D<T, Descriptor>(
                        surfaceTension, incompressibleModel),
                    group.getBoundingBox(), freeSurfaceArgs, pl);
            } else {
                std::vector<MultiBlock3D *> args(freeSurfaceArgs);
                args.push_back(surfaceTensionField);

                integrateProcessingFunctional(
                    new FreeSurfaceAddSurfaceTensionFromScalarField3D<T, Descriptor>(
                        incompressibleModel),
                    group.getBoundingBox(), args, pl);
            }
        }

        integrateProcessingFunctional(
            new FreeSurfaceStabilize3D<T, Descriptor>(), group.getBoundingBox(), freeSurfaceArgs,
            pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), group.getBoundingBox(),
            freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(
                rhoDefault, emptyNodeDynamics->clone()),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(dynamics->clone(), force),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(
                rhoDefault, emptyNodeDynamics->clone()),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceEqualMassExcessReDistributionAndComputationOfLostMass3D<T, Descriptor>(
                reductionData),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        {
            std::vector<MultiBlock3D *> args;
            args.push_back(&group.get(fsPrefix + "flag"));
            args.push_back(&group.get(fsPrefix + "mass"));

            integrateProcessingFunctional(
                new FreeSurfaceComputeReductionsPerProcess3D<T>(reductionData),
                group.getBoundingBox(), fluidLattice, args, pl);
        }

        /***** New level ******/
        pl++;

        {
            std::vector<MultiBlock3D *> args;
            args.push_back(&fluidLattice);

#ifdef PLB_MPI_PARALLEL
            createReductionCommunicator(fluidLattice);
            integrateProcessingFunctional(
                new FreeSurfaceComputeReductions3D(reductionData, reductionCommunicator),
                group.getBoundingBox(), args, pl);
#else
            integrateProcessingFunctional(
                new FreeSurfaceComputeSerialReductions3D(reductionData), group.getBoundingBox(),
                args, pl);
#endif
        }

        /***** New level ******/
        pl++;

        {
            std::vector<MultiBlock3D *> args;
            args.push_back(&fluidLattice);

            integrateProcessingFunctional(
                new FreeSurfaceResetReductionData3D(reductionData), group.getBoundingBox(), args,
                pl);
        }

        /***** New level ******/
        pl++;

        if (useForce) {
            integrateProcessingFunctional(
                new FreeSurfaceAddExternalForce3D<T, Descriptor>(rhoDefault),
                group.getBoundingBox(), freeSurfaceArgs, pl);
        } else if (forceProcessor) {
            integrateProcessingFunctional(
                forceProcessor, group.getBoundingBox(), fluidLattice, blocks, pl);
        }

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceLostMassReDistribution3D<T, Descriptor>(reductionData),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        if (collideAndStreamAtTheEndOfCycle) {
            pl++;

            integrateProcessingFunctional(
                new ExternalRhoJcollideAndStream3D<T, Descriptor>, group.getBoundingBox(),
                rhoBarJparam, pl);
        }
    }

    void createFreeSurfaceProcessors(
        plint initialProcessorLevel, bool collideAndStreamAtTheEndOfCycle,
        BoxProcessingFunctional3D *forceProcessor, std::vector<MultiBlock3D *> blocks)
    {
        MultiScalarField3D<T> *surfaceTensionField = 0;
        createFreeSurfaceProcessors(
            initialProcessorLevel, collideAndStreamAtTheEndOfCycle, forceProcessor, blocks,
            surfaceTensionField);
    }

    void createFreeSurfaceProcessors(
        plint initialProcessorLevel = 0, bool collideAndStreamAtTheEndOfCycle = false,
        MultiScalarField3D<T> *surfaceTensionField = 0)
    {
        BoxProcessingFunctional3D *forceProcessor = 0;
        std::vector<MultiBlock3D *> blocks;
        createFreeSurfaceProcessors(
            initialProcessorLevel, collideAndStreamAtTheEndOfCycle, forceProcessor, blocks,
            surfaceTensionField);
    }

    template <class VelFunction>
    void createFreeSurfaceProcessorsForImmersedWalls(
        plint numIBIterations, std::vector<Array<T, 3> > const &vertices,
        std::vector<T> const &areas, std::vector<int> const &flags, VelFunction velFunction,
        int repelInterface, plint initialProcessorLevel, bool collideAndStreamAtTheEndOfCycle,
        BoxProcessingFunctional3D *forceProcessor, std::vector<MultiBlock3D *> blocks,
        MultiScalarField3D<T> *surfaceTensionField)
    {
        PLB_ASSERT(initialProcessorLevel >= 0);

        bool useForce = !util::isZero(norm(force));
        PLB_ASSERT(!(useForce && forceProcessor));

        bool useSurfaceTension = (surfaceTensionField != 0 || !util::isZero(surfaceTension));

        initializeInterfaceLists3D<T, Descriptor>(group.getContainer(fsPrefix + "helperLists"));
        // setToConstant(group.getScalar<int>("flag"), group.getBoundingBox(),
        // (int)freeSurfaceFlag::empty); setToConstant(group.getScalar<T>("outsideDensity"),
        // group.getBoundingBox(), rhoDefault);

        MultiBlockLattice3D<T, Descriptor> &fluidLattice =
            (lattice ? *lattice : group.getLattice<T, Descriptor>(fluidname));

        plint pl = initialProcessorLevel;  // Processor level.

        /***** Initial level ******/
        if (pl == 0 && !collideAndStreamAtTheEndOfCycle) {
            integrateProcessingFunctional(
                new ExternalRhoJcollideAndStream3D<T, Descriptor>, group.getBoundingBox(),
                rhoBarJparam, pl);
        }

        integrateProcessingFunctional(
            new FreeSurfaceComputeNormals3D<T, Descriptor>, group.getBoundingBox(), freeSurfaceArgs,
            pl);

        /***** New level ******/
        pl++;

        if (useSurfaceTension) {
            if (contactAngleFunction == 0) {
                integrateProcessingFunctional(
                    new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngle),
                    group.getBoundingBox(), freeSurfaceArgs, pl);
            } else {
                integrateProcessingFunctional(
                    new FreeSurfaceComputeCurvature3D<T, Descriptor>(contactAngleFunction),
                    group.getBoundingBox(), freeSurfaceArgs, pl);
            }

            // To change to the curvature calculation with height functions, uncomment the next data
            // processor and comment out the two previous ones. If only the next data processor is
            // used and there is no surface tension, the normals are not computed at all. Be careful
            // if you intent to use the normals and do not have the surface tension algorithm
            // enabled.
            // integrateProcessingFunctional (
            //        new FreeSurfaceGeometry3D<T,Descriptor>(contactAngle),
            //        group.getBoundingBox(), freeSurfaceArgs, pl );
        }

        integrateProcessingFunctional(
            new FreeSurfaceMassChange3D<T, Descriptor>, group.getBoundingBox(), freeSurfaceArgs,
            pl);

        integrateProcessingFunctional(
            new FreeSurfaceCompletion3D<T, Descriptor>, group.getBoundingBox(), freeSurfaceArgs,
            pl);

        integrateProcessingFunctional(
            new FreeSurfaceMacroscopicWithoutLostMassReDistribution3D<T, Descriptor>(
                incompressibleModel),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        if (useSurfaceTension) {
            if (surfaceTensionField == 0) {
                integrateProcessingFunctional(
                    new FreeSurfaceAddSurfaceTension3D<T, Descriptor>(
                        surfaceTension, incompressibleModel),
                    group.getBoundingBox(), freeSurfaceArgs, pl);
            } else {
                std::vector<MultiBlock3D *> args(freeSurfaceArgs);
                args.push_back(surfaceTensionField);

                integrateProcessingFunctional(
                    new FreeSurfaceAddSurfaceTensionFromScalarField3D<T, Descriptor>(
                        incompressibleModel),
                    group.getBoundingBox(), args, pl);
            }
        }

        integrateProcessingFunctional(
            new FreeSurfaceStabilize3D<T, Descriptor>(), group.getBoundingBox(), freeSurfaceArgs,
            pl);

        std::vector<MultiBlock3D *> immersedWallDataArgs;
        immersedWallDataArgs.push_back(&group.get(fsPrefix + "ibm_container"));
        integrateProcessingFunctional(
            new InstantiateImmersedWallDataWithIndexedTagging3D<T>(vertices, areas, flags),
            group.getBoundingBox(), fluidLattice, immersedWallDataArgs, pl);

        /***** New level ******/
        if (repelInterface == 2) {
            pl++;

            std::vector<MultiBlock3D *> tmpProtectionArgs;
            tmpProtectionArgs.push_back(&group.get(fsPrefix + "flag"));
            tmpProtectionArgs.push_back(&group.get(fsPrefix + "ibm_container"));
            integrateProcessingFunctional(
                new TemporarilyProtectImmersedWalls3D<T>(), group.getBoundingBox(), fluidLattice,
                tmpProtectionArgs, pl);
        }

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceComputeInterfaceLists3D<T, Descriptor>(), group.getBoundingBox(),
            freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>(
                rhoDefault, emptyNodeDynamics->clone()),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        integrateProcessingFunctional(
            new FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>(dynamics->clone(), force),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>(
                rhoDefault, emptyNodeDynamics->clone()),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceEqualMassExcessReDistributionAndComputationOfLostMass3D<T, Descriptor>(
                reductionData),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        pl++;

        {
            std::vector<MultiBlock3D *> args;
            args.push_back(&group.get(fsPrefix + "flag"));
            args.push_back(&group.get(fsPrefix + "mass"));

            integrateProcessingFunctional(
                new FreeSurfaceComputeReductionsPerProcess3D<T>(reductionData),
                group.getBoundingBox(), fluidLattice, args, pl);
        }

        /***** New level ******/
        pl++;

        {
            std::vector<MultiBlock3D *> args;
            args.push_back(&fluidLattice);

#ifdef PLB_MPI_PARALLEL
            createReductionCommunicator(fluidLattice);
            integrateProcessingFunctional(
                new FreeSurfaceComputeReductions3D(reductionData, reductionCommunicator),
                group.getBoundingBox(), args, pl);
#else
            integrateProcessingFunctional(
                new FreeSurfaceComputeSerialReductions3D(reductionData), group.getBoundingBox(),
                args, pl);
#endif
        }

        /***** New level ******/
        pl++;

        {
            std::vector<MultiBlock3D *> args;
            args.push_back(&fluidLattice);

            integrateProcessingFunctional(
                new FreeSurfaceResetReductionData3D(reductionData), group.getBoundingBox(), args,
                pl);
        }

        /***** New level ******/

        for (int i = 0; i < numIBIterations; i++) {
            pl++;
            T tau = (T)1 / dynamics->getOmega();
            std::vector<MultiBlock3D *> args;
            args.push_back(&group.get(rhoBarJprefix + "rhoBar"));
            args.push_back(&group.get(rhoBarJprefix + "j"));
            args.push_back(&group.get(fsPrefix + "ibm_container"));
            integrateProcessingFunctional(
                new IndexedInamuroIteration3D<T, VelFunction>(
                    velFunction, tau, incompressibleModel),
                group.getBoundingBox(), fluidLattice, args, pl);
        }

        if (repelInterface == 1 || repelInterface == 3) {
            bool strongRepelling = repelInterface == 1 ? false : true;

            // TODO: This data processor changes the momentum in the vicinity of the
            //       immersed boundary. Maybe it affects strongly the measurement of
            //       force and torque on the immersed surface. Needs to be checked.
            pl++;
            std::vector<MultiBlock3D *> args;
            args.push_back(&group.get(rhoBarJprefix + "rhoBar"));
            args.push_back(&group.get(rhoBarJprefix + "j"));
            args.push_back(&group.get(fsPrefix + "flag"));
            args.push_back(&group.get(fsPrefix + "ibm_container"));
            integrateProcessingFunctional(
                new RepelInterfaceFromImmersedWalls3D<T, VelFunction>(
                    velFunction, rhoDefault, strongRepelling),
                group.getBoundingBox(), fluidLattice, args, pl);
        }

        /***** New level ******/
        if (repelInterface == 2) {
            pl++;

            std::vector<MultiBlock3D *> rmProtectionArgs;
            rmProtectionArgs.push_back(&group.get(fsPrefix + "flag"));
            integrateProcessingFunctional(
                new RemoveProtectionFromImmersedWalls3D<T>(), group.getBoundingBox(), fluidLattice,
                rmProtectionArgs, pl);
        }

        /***** New level ******/
        pl++;

        if (useForce) {
            integrateProcessingFunctional(
                new FreeSurfaceAddExternalForce3D<T, Descriptor>(rhoDefault),
                group.getBoundingBox(), freeSurfaceArgs, pl);
        } else if (forceProcessor) {
            integrateProcessingFunctional(
                forceProcessor, group.getBoundingBox(), fluidLattice, blocks, pl);
        }

        /***** New level ******/
        pl++;

        integrateProcessingFunctional(
            new FreeSurfaceLostMassReDistribution3D<T, Descriptor>(reductionData),
            group.getBoundingBox(), freeSurfaceArgs, pl);

        /***** New level ******/
        if (collideAndStreamAtTheEndOfCycle) {
            pl++;

            integrateProcessingFunctional(
                new ExternalRhoJcollideAndStream3D<T, Descriptor>, group.getBoundingBox(),
                rhoBarJparam, pl);
        }
    }

    template <class VelFunction>
    void createFreeSurfaceProcessorsForImmersedWalls(
        plint numIBIterations, std::vector<Array<T, 3> > const &vertices,
        std::vector<T> const &areas, std::vector<int> const &flags, VelFunction velFunction,
        int repelInterface, plint initialProcessorLevel = 0,
        bool collideAndStreamAtTheEndOfCycle = false,
        MultiScalarField3D<T> *surfaceTensionField = 0)
    {
        BoxProcessingFunctional3D *forceProcessor = 0;
        std::vector<MultiBlock3D *> blocks;

        createFreeSurfaceProcessorsForImmersedWalls<VelFunction>(
            numIBIterations, vertices, areas, flags, velFunction, repelInterface,
            initialProcessorLevel, collideAndStreamAtTheEndOfCycle, forceProcessor, blocks,
            surfaceTensionField);
    }

    void defaultInitialize(
        bool useConstRho = true, bool useZeroMomentum = true, bool initializeCell = true)
    {
        applyProcessingFunctional(
            new DefaultInitializeFreeSurface3D<T, Descriptor>(
                dynamics->clone(), emptyNodeDynamics->clone(), force, rhoDefault, useConstRho,
                useZeroMomentum, initializeCell),
            group.getBoundingBox(), freeSurfaceArgs);
    }

    void partiallyDefaultInitialize(
        bool useConstRho = true, bool useZeroMomentum = true, bool initializeCell = true)
    {
        applyProcessingFunctional(
            new PartiallyDefaultInitializeFreeSurface3D<T, Descriptor>(
                dynamics->clone(), emptyNodeDynamics->clone(), force, rhoDefault, useConstRho,
                useZeroMomentum, initializeCell),
            group.getBoundingBox(), freeSurfaceArgs);
    }

    Group3D &getGroup()
    {
        return group;
    }

    Group3D const &getGroup() const
    {
        return group;
    }

    MultiBlockLattice3D<T, Descriptor> &getLattice()
    {
        MultiBlockLattice3D<T, Descriptor> &fluidLattice =
            (lattice ? *lattice : group.getLattice<T, Descriptor>(fluidname));
        return fluidLattice;
    }

    MultiBlockLattice3D<T, Descriptor> const &getLattice() const
    {
        MultiBlockLattice3D<T, Descriptor> &fluidLattice =
            (lattice ? *lattice : group.getLattice<T, Descriptor>(fluidname));
        return fluidLattice;
    }

    Dynamics<T, Descriptor> const &getDynamics() const
    {
        return *dynamics;
    }

    MultiScalarField3D<int> &getFlag()
    {
        return group.getScalar<int>(fsPrefix + "flag");
    }

    MultiScalarField3D<int> const &getFlag() const
    {
        return group.getScalar<int>(fsPrefix + "flag");
    }

    MultiScalarField3D<T> &getMass()
    {
        return group.getScalar<T>(fsPrefix + "mass");
    }

    MultiScalarField3D<T> const &getMass() const
    {
        return group.getScalar<T>(fsPrefix + "mass");
    }

    MultiScalarField3D<T> &getVolumeFraction()
    {
        return group.getScalar<T>(fsPrefix + "volumeFraction");
    }

    MultiScalarField3D<T> const &getVolumeFraction() const
    {
        return group.getScalar<T>(fsPrefix + "volumeFraction");
    }

    MultiScalarField3D<T> &getRhoBar()
    {
        return group.getScalar<T>(rhoBarJprefix + "rhoBar");
    }

    MultiScalarField3D<T> const &getRhoBar() const
    {
        return group.getScalar<T>(rhoBarJprefix + "rhoBar");
    }

    MultiTensorField3D<T, 3> &getJ()
    {
        return group.getTensor<T, 3>(rhoBarJprefix + "j");
    }

    MultiTensorField3D<T, 3> const &getJ() const
    {
        return group.getTensor<T, 3>(rhoBarJprefix + "j");
    }

    T getRhoDefault() const
    {
        return rhoDefault;
    }

    T getSurfaceTension() const
    {
        return surfaceTension;
    }

    T getContactAngle() const
    {
        return contactAngle;
    }

    ContactAngleFunction getContactAngleFunction() const
    {
        return contactAngleFunction;
    }

    Array<T, Descriptor<T>::ExternalField::sizeOfForce> const &getForce() const
    {
        return force;
    }

    bool getIncompressibleModel() const
    {
        return incompressibleModel;
    }

    FreeSurfaceReductionData &getReductionData()
    {
        return reductionData;
    }

    FreeSurfaceReductionData const &getReductionData() const
    {
        return reductionData;
    }

#ifdef PLB_MPI_PARALLEL
    MPI_Comm &getReductionCommunicator()
    {
        return reductionCommunicator;
    }

    MPI_Comm const &getReductionCommunicator() const
    {
        return reductionCommunicator;
    }
#endif

    T getLostMass() const
    {
        return reductionData.lostMass;
    }

    T getTotalMass() const
    {
        return reductionData.totalMass;
    }

    plint getNumInterfaceCells() const
    {
        return reductionData.numInterfaceCells;
    }

private:
    // TODO: maybe replace ASSERTS with exceptions.
    bool fieldExists(std::string name, [[maybe_unused]] plint envelopeWidth)
    {
        if (group.hasBlock(name)) {
            PLB_ASSERT(
                group.get(name).getMultiBlockManagement().getEnvelopeWidth() >= envelopeWidth);
            return true;
        }
        return false;
    }

#ifdef PLB_MPI_PARALLEL
    void createReductionCommunicator(MultiBlock3D &multiBlock)
    {
        SparseBlockStructure3D const &sparseBlockStructure =
            multiBlock.getMultiBlockManagement().getSparseBlockStructure();
        ThreadAttribution const &threadAttribution =
            multiBlock.getMultiBlockManagement().getThreadAttribution();

        std::map<plint, Box3D> const &bulks = sparseBlockStructure.getBulks();
        std::set<int> processIds;
        std::map<plint, Box3D>::const_iterator bulkIt = bulks.begin();
        for (; bulkIt != bulks.end(); ++bulkIt) {
            plint blockId = bulkIt->first;
            (void)processIds.insert(threadAttribution.getMpiProcess(blockId));
        }

        std::vector<int> processIdsVec;
        std::set<int>::const_iterator processIdsIt = processIds.begin();
        for (; processIdsIt != processIds.end(); ++processIdsIt) {
            processIdsVec.push_back(*processIdsIt);
        }
        PLB_ASSERT(!processIdsVec.empty());

        MPI_Group globalGroup;
        (void)MPI_Comm_group(global::mpi().getGlobalCommunicator(), &globalGroup);
        MPI_Group reductionGroup;
        (void)MPI_Group_incl(globalGroup, processIdsVec.size(), &processIdsVec[0], &reductionGroup);

        (void)MPI_Comm_create(
            global::mpi().getGlobalCommunicator(), reductionGroup, &reductionCommunicator);
    }
#endif
private:
    Group3D &group;
    MultiBlockLattice3D<T, Descriptor>
        *lattice;  // Pointer to a fluid lattice if not contained in the group.
    Dynamics<T, Descriptor> *dynamics;
    Dynamics<T, Descriptor> *emptyNodeDynamics;
    bool incompressibleModel;
    T rhoDefault, surfaceTension, contactAngle;
    ContactAngleFunction contactAngleFunction;
    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force;
    std::string fluidname, fsPrefix, rhoBarJprefix;
    bool hasImmersedWalls;
    FreeSurfaceReductionData reductionData;
#ifdef PLB_MPI_PARALLEL
    MPI_Comm reductionCommunicator;
#endif
public:
    std::vector<MultiBlock3D *> freeSurfaceArgs;
    std::vector<MultiBlock3D *> rhoBarJparam;
};

}  // namespace plb

#endif  // FREE_SURFACE_MODEL_3D_H
