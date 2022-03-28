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

#ifndef FREE_SURFACE_UTIL_3D_H
#define FREE_SURFACE_UTIL_3D_H

#include <set>
#include <string>
#include <vector>

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "finiteDifference/fdWrapper3D.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "latticeBoltzmann/indexTemplates.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

/// Constants used in a free surface flag matrix for cell tagging.
namespace freeSurfaceFlag {
enum Flag {
    empty = 0,
    interface = 1,
    fluid = 2,
    wall = 4,
    protect = 5,
    temporarilyProtect = 6,
    protectEmpty = 7,
    slipWall = 8
};
inline std::string flagToString(int flag)
{
    switch (flag) {
    case empty:
        return "empty";
    case interface:
        return "interface";
    case fluid:
        return "fluid";
    case wall:
        return "wall";
    case protect:
        return "protect";
    case temporarilyProtect:
        return "temporarilyProtect";
    case protectEmpty:
        return "protectEmpty";
    case slipWall:
        return "slipWall";
    default:
        PLB_ASSERT(false);
    }
    return std::string();
}
inline Flag invert(int flag)
{
    switch (flag) {
    case empty:
        return fluid;
    case interface:
        return interface;
    case fluid:
        return empty;
    case wall:
        return wall;
    case protect:
        return protect;
    case temporarilyProtect:
        return temporarilyProtect;
    case protectEmpty:
        return protectEmpty;
    case slipWall:
        return slipWall;
    default:
        PLB_ASSERT(false);
    }
    return (Flag)(-1);
}
inline bool isWet(int flag)
{
    return flag == interface || flag == fluid || flag == protect || flag == temporarilyProtect;
}
inline bool isFullWet(int flag)
{
    return flag == fluid || flag == protect || flag == temporarilyProtect;
}
inline bool isProtected(int flag)
{
    return flag == protect || flag == temporarilyProtect;
}
inline bool isEmpty(int flag)
{
    return flag == empty || flag == protectEmpty;
}
inline bool isAnyWall(int flag)
{
    return flag == wall || flag == slipWall;
}
inline bool isWall(int flag)
{
    return flag == wall;
}
inline bool isSlipWall(int flag)
{
    return flag == slipWall;
}
}  // namespace freeSurfaceFlag

/// Create a parameter-list for most free-surface data processors.
template <typename T, template <typename U> class Descriptor>
std::vector<MultiBlock3D *> aggregateFreeSurfaceParams(
    MultiBlockLattice3D<T, Descriptor> &fluid, MultiScalarField3D<T> &rhoBar,
    MultiTensorField3D<T, 3> &j, MultiScalarField3D<T> &mass, MultiScalarField3D<T> &volumeFraction,
    MultiScalarField3D<int> &flag, MultiTensorField3D<T, 3> &normal,
    MultiContainerBlock3D &interfaceLists, MultiScalarField3D<T> &curvature,
    MultiScalarField3D<T> &outsideDensity)
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

    return aggregation;
}

/// Data structure for holding lists of cells along the free surface in an AtomicContainerBlock.
template <typename T, template <typename U> class Descriptor>
struct InterfaceLists : public ContainerBlockData {
    typedef Array<plint, Descriptor<T>::d> Node;
    /// Holds all nodes which have excess mass.
    std::map<Node, T> massExcess;
    /// Holds all nodes that need to change status from interface to fluid.
    std::set<Node> interfaceToFluid;
    /// Holds all nodes that need to change status from interface to empty.
    std::set<Node> interfaceToEmpty;
    /// Holds all nodes that need to change status from empty to interface.
    std::set<Node> emptyToInterface;

    virtual InterfaceLists<T, Descriptor> *clone() const
    {
        return new InterfaceLists<T, Descriptor>(*this);
    }
};

// Reductions are performed always in double precision.
struct FreeSurfaceReductionData {
    FreeSurfaceReductionData() :
        localLostMass((double)0),
        lostMass((double)0),
        localTotalMass((double)0),
        totalMass((double)0),
        localNumInterfaceCells((plint)0),
        numInterfaceCells((plint)0),
        used(false)
    { }

    void reset()
    {
        localLostMass = (double)0;
        localTotalMass = (double)0;
        localNumInterfaceCells = (plint)0;
        used = false;
    }

    double localLostMass, lostMass;
    double localTotalMass, totalMass;
    plint localNumInterfaceCells, numInterfaceCells;
    bool used;
};

/// A wrapper offering convenient access to the free-surface data provided to
/// data processors. Avoids verbous casting, asserting, etc.
template <typename T, template <typename U> class Descriptor>
class FreeSurfaceProcessorParam3D {
public:
    typedef typename InterfaceLists<T, Descriptor>::Node Node;
    FreeSurfaceProcessorParam3D(std::vector<AtomicBlock3D *> &atomicBlocks)
    {
        PLB_ASSERT(atomicBlocks.size() >= 10);

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

        interfaceLists_ =
            dynamic_cast<InterfaceLists<T, Descriptor> *>(containerInterfaceLists_->getData());
        // PLB_ASSERT(interfaceLists_); // This assertion must be commented out to work with
        // TwoPhase.

        curvature_ = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[8]);
        PLB_ASSERT(curvature_);

        outsideDensity_ = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[9]);
        PLB_ASSERT(outsideDensity_);

        absoluteOffset = fluid_->getLocation();
        relativeOffsetRhoBar = computeRelativeDisplacement(*fluid_, *rhoBar_);
        relativeOffsetJ = computeRelativeDisplacement(*fluid_, *j_);
        relativeOffsetMass = computeRelativeDisplacement(*fluid_, *mass_);
        relativeOffsetVF = computeRelativeDisplacement(*fluid_, *volumeFraction_);
        relativeOffsetFS = computeRelativeDisplacement(*fluid_, *flag_);
        relativeOffsetNormal = computeRelativeDisplacement(*fluid_, *normal_);
        relativeOffsetC = computeRelativeDisplacement(*fluid_, *curvature_);
        relativeOffsetOD = computeRelativeDisplacement(*fluid_, *outsideDensity_);
    }

    Cell<T, Descriptor> &cell(plint iX, plint iY, plint iZ)
    {
        return fluid_->get(iX, iY, iZ);
    }
    T &mass(plint iX, plint iY, plint iZ)
    {
        return mass_->get(
            iX + relativeOffsetMass.x, iY + relativeOffsetMass.y, iZ + relativeOffsetMass.z);
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
    void setMomentum(plint iX, plint iY, plint iZ, Array<T, 3> const &momentum)
    {
        j_->get(iX + relativeOffsetJ.x, iY + relativeOffsetJ.y, iZ + relativeOffsetJ.z) = momentum;
    }
    Array<T, 3> getMomentum(plint iX, plint iY, plint iZ)
    {
        return j_->get(iX + relativeOffsetJ.x, iY + relativeOffsetJ.y, iZ + relativeOffsetJ.z);
    }
    T getDensity(plint iX, plint iY, plint iZ)
    {
        return Descriptor<T>::fullRho(rhoBar_->get(
            iX + relativeOffsetRhoBar.x, iY + relativeOffsetRhoBar.y, iZ + relativeOffsetRhoBar.z));
    }
    void setDensity(plint iX, plint iY, plint iZ, T rho)
    {
        rhoBar_->get(
            iX + relativeOffsetRhoBar.x, iY + relativeOffsetRhoBar.y, iZ + relativeOffsetRhoBar.z) =
            Descriptor<T>::rhoBar(rho);
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
    plint getNumInterfaceCells() const
    {
        return fluid_->getInternalStatistics().getIntSum(0);
    }

    T smooth(ScalarField3D<T> const &scalar, Dot3D const &ofs, plint iX, plint iY, plint iZ)
    {
        using namespace freeSurfaceFlag;

        if (isAnyWall(flag_->get(
                iX + relativeOffsetFS.x, iY + relativeOffsetFS.y, iZ + relativeOffsetFS.z)))
        {
            return (scalar.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z));
        }

        T val = 0.0;
        int n = 0;
        for (int i = -1; i < 2; i++) {
            plint nextX = iX + i;
            for (int j = -1; j < 2; j++) {
                plint nextY = iY + j;
                for (int k = -1; k < 2; k++) {
                    plint nextZ = iZ + k;
                    if (!(i == 0 && j == 0 && k == 0)
                        && !isAnyWall(flag_->get(
                            nextX + relativeOffsetFS.x, nextY + relativeOffsetFS.y,
                            nextZ + relativeOffsetFS.z)))
                    {
                        n++;
                        val += scalar.get(nextX + ofs.x, nextY + ofs.y, nextZ + ofs.z);
                    }
                }
            }
        }
        if (n != 0) {
            val /= (T)n;
        } else {
            val = scalar.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z);
        }

        return (val);
    }

    // The following function can be called with a lattice descriptor "Descriptor2", which is
    // different from the lattice descriptor "Descriptor" used in the simulation.
    template <typename T2, template <typename U2> class Descriptor2>
    T lbmSmooth(ScalarField3D<T> const &scalar, Dot3D const &ofs, plint iX, plint iY, plint iZ)
    {
        typedef Descriptor2<T2> D;
        using namespace freeSurfaceFlag;

        if (isAnyWall(flag_->get(
                iX + relativeOffsetFS.x, iY + relativeOffsetFS.y, iZ + relativeOffsetFS.z)))
        {
            return (scalar.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z));
        }

        T val = 0.0, sum = 0.0;
        for (plint iPop = 1; iPop < D::q; iPop++) {
            plint nextX = iX + D::c[iPop][0];
            plint nextY = iY + D::c[iPop][1];
            plint nextZ = iZ + D::c[iPop][2];
            sum += D::t[iPop];

            T nextVal = 0.0;

            // First, extrapolate the scalar field on the wall (if necessary).
            if (isAnyWall(flag_->get(
                    nextX + relativeOffsetFS.x, nextY + relativeOffsetFS.y,
                    nextZ + relativeOffsetFS.z)))
            {
                T locVal = scalar.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z);
                plint opp = indexTemplates::opposite<D>(iPop);
                plint prevX = iX + D::c[opp][0];
                plint prevY = iY + D::c[opp][1];
                plint prevZ = iZ + D::c[opp][2];
                if (isAnyWall(flag_->get(
                        prevX + relativeOffsetFS.x, prevY + relativeOffsetFS.y,
                        prevZ + relativeOffsetFS.z)))
                {
                    nextVal = locVal;
                } else {
                    T prevVal = scalar.get(prevX + ofs.x, prevY + ofs.y, prevZ + ofs.z);
                    nextVal = (T)2 * locVal - prevVal;
                }
            } else {
                nextVal = scalar.get(nextX + ofs.x, nextY + ofs.y, nextZ + ofs.z);
            }

            val += D::t[iPop] * nextVal;
        }
        val /= sum;

        return (val);
    }

    // The following function can be called with a lattice descriptor "Descriptor2", which is
    // different from the lattice descriptor "Descriptor" used in the simulation.
    template <typename T2, template <typename U2> class Descriptor2>
    Array<T, 3> lbmSmooth(
        TensorField3D<T, 3> const &tensor, Dot3D const &ofs, plint iX, plint iY, plint iZ)
    {
        typedef Descriptor2<T2> D;
        using namespace freeSurfaceFlag;

        if (isAnyWall(flag_->get(
                iX + relativeOffsetFS.x, iY + relativeOffsetFS.y, iZ + relativeOffsetFS.z)))
        {
            return (tensor.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z));
        }

        T sum = 0.0;
        Array<T, 3> val((T)0, (T)0, (T)0);
        for (plint iPop = 1; iPop < D::q; iPop++) {
            plint nextX = iX + D::c[iPop][0];
            plint nextY = iY + D::c[iPop][1];
            plint nextZ = iZ + D::c[iPop][2];
            sum += D::t[iPop];

            Array<T, 3> nextVal((T)0, (T)0, (T)0);

            // First, extrapolate the tensor field on the wall (if necessary).
            if (isAnyWall(flag_->get(
                    nextX + relativeOffsetFS.x, nextY + relativeOffsetFS.y,
                    nextZ + relativeOffsetFS.z)))
            {
                Array<T, 3> locVal = tensor.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z);
                plint opp = indexTemplates::opposite<D>(iPop);
                plint prevX = iX + D::c[opp][0];
                plint prevY = iY + D::c[opp][1];
                plint prevZ = iZ + D::c[opp][2];
                if (isAnyWall(flag_->get(
                        prevX + relativeOffsetFS.x, prevY + relativeOffsetFS.y,
                        prevZ + relativeOffsetFS.z)))
                {
                    nextVal = locVal;
                } else {
                    Array<T, 3> prevVal = tensor.get(prevX + ofs.x, prevY + ofs.y, prevZ + ofs.z);
                    nextVal = (T)2 * locVal - prevVal;
                }
            } else {
                nextVal = tensor.get(nextX + ofs.x, nextY + ofs.y, nextZ + ofs.z);
            }

            val += D::t[iPop] * nextVal;
        }
        val /= sum;

        return (val);
    }

    // In the free-surface algorithm there is the need to compute derivatives with finite
    // differences. The following function, computes the x/y/z widths and positions for finite
    // difference stencils depending on the local point in the simulation domain. The stencils are
    // tuned so that lattice nodes inside walls are excluded. The maximum stencil width is always an
    // odd number of the form: 2 * h + 1, with h being less or equal to the envelope width of the
    // flag matrix.
    void getFdStencilWidthsAndPositions(
        plint iX, plint iY, plint iZ, plint h, Array<int, 3> &widths, Array<int, 3> &positions)
    {
        using namespace freeSurfaceFlag;

        plint x = iX + relativeOffsetFS.x;
        plint y = iY + relativeOffsetFS.y;
        plint z = iZ + relativeOffsetFS.z;

        // PLB_ASSERT(!isAnyWall(flag_->get(x, y, z)));

        int left, right;

        // x-direction.

        left = 0;
        for (plint d = 1; d <= h; d++, left++) {
            if (isAnyWall(flag_->get(x - d, y, z))) {
                break;
            }
        }
        right = 0;
        for (plint d = 1; d <= h; d++, right++) {
            if (isAnyWall(flag_->get(x + d, y, z))) {
                break;
            }
        }

        widths[0] = left + right + 1;
        positions[0] = left;

        // y-direction.

        left = 0;
        for (plint d = 1; d <= h; d++, left++) {
            if (isAnyWall(flag_->get(x, y - d, z))) {
                break;
            }
        }
        right = 0;
        for (plint d = 1; d <= h; d++, right++) {
            if (isAnyWall(flag_->get(x, y + d, z))) {
                break;
            }
        }

        widths[1] = left + right + 1;
        positions[1] = left;

        // z-direction.

        left = 0;
        for (plint d = 1; d <= h; d++, left++) {
            if (isAnyWall(flag_->get(x, y, z - d))) {
                break;
            }
        }
        right = 0;
        for (plint d = 1; d <= h; d++, right++) {
            if (isAnyWall(flag_->get(x, y, z + d))) {
                break;
            }
        }

        widths[2] = left + right + 1;
        positions[2] = left;
    }

    // Compute the gradient of a scalar field with finite differences excluding the wall cells.
    Array<T, 3> computeGradient(
        ScalarField3D<T> const &scalar, Dot3D const &ofs, plint h, plint iX, plint iY, plint iZ)
    {
        using namespace freeSurfaceFlag;

        // PLB_ASSERT(h >= 0 && h <= 3);
        // PLB_ASSERT(!isAnyWall(flag_->get(iX + relativeOffsetFS.x, iY + relativeOffsetFS.y, iZ +
        // relativeOffsetFS.z)));

        Array<int, 3> widths, positions;
        getFdStencilWidthsAndPositions(iX, iY, iZ, h, widths, positions);

        plint i = iX + ofs.x;
        plint j = iY + ofs.y;
        plint k = iZ + ofs.z;

        Array<T, 3> gradient;

        // Template parameters:
        // 1: because we compute first order derivatives
        // 7: this is at least (2 * h + 1), with "h" being at most equal to the envelope width of
        // the flag matrix
        //    We set this to its current possible maximum value 7, to be more efficient in case we
        //    want to use high-order finite differences. This must be changed if "h" is ever
        //    expected to be larger than 3.
        gradient[0] = computeScalarXderivative<T, 1, 7>(scalar, widths[0], positions[0], i, j, k);
        gradient[1] = computeScalarYderivative<T, 1, 7>(scalar, widths[1], positions[1], i, j, k);
        gradient[2] = computeScalarZderivative<T, 1, 7>(scalar, widths[2], positions[2], i, j, k);

        return (gradient);
    }

    // Compute the divergence of a tensor field with finite differences excluding the wall cells.
    T computeDivergence(
        TensorField3D<T, 3> const &tensor, Dot3D const &ofs, plint h, plint iX, plint iY, plint iZ)
    {
        using namespace freeSurfaceFlag;

        // PLB_ASSERT(h >= 0 && h <= 3);
        // PLB_ASSERT(!isAnyWall(flag_->get(iX + relativeOffsetFS.x, iY + relativeOffsetFS.y, iZ +
        // relativeOffsetFS.z)));

        Array<int, 3> widths, positions;
        getFdStencilWidthsAndPositions(iX, iY, iZ, h, widths, positions);

        plint i = iX + ofs.x;
        plint j = iY + ofs.y;
        plint k = iZ + ofs.z;

        // Template parameters:
        // 3: because "tensor" contains 3D vectors
        // 1: because we compute first order derivatives
        // 7: this is at least (2 * h + 1), with "h" being at most equal to the envelope width of
        // the flag matrix
        //    We set this to its current possible maximum value 7, to be more efficient in case we
        //    want to use high-order finite differences. This must be changed if "h" is ever
        //    expected to be larger than 3.
        Array<T, 3> d_dx =
            computeTensorXderivative<T, 3, 1, 7>(tensor, widths[0], positions[0], i, j, k);
        Array<T, 3> d_dy =
            computeTensorYderivative<T, 3, 1, 7>(tensor, widths[1], positions[1], i, j, k);
        Array<T, 3> d_dz =
            computeTensorZderivative<T, 3, 1, 7>(tensor, widths[2], positions[2], i, j, k);

        T divergence = d_dx[0] + d_dy[1] + d_dz[2];

        return (divergence);
    }

    // Compute the gradient of a scalar field "the lattice Boltzmann way".
    // With this method the scalar field values are first extrapolated on the neighboring wall cells
    // (if necessary). So, in contrast to the above finite difference implementation, the wall cells
    // are not excluded here. The following function can be called with a lattice descriptor
    // "Descriptor2", which is different from the lattice descriptor "Descriptor" used in the
    // simulation.
    template <typename T2, template <typename U2> class Descriptor2>
    Array<T, 3> lbmComputeGradient(
        ScalarField3D<T> const &scalar, Dot3D const &ofs, plint iX, plint iY, plint iZ)
    {
        typedef Descriptor2<T2> D;
        using namespace freeSurfaceFlag;

        // PLB_ASSERT(!isAnyWall(flag_->get(iX + relativeOffsetFS.x, iY + relativeOffsetFS.y, iZ +
        // relativeOffsetFS.z)));

        Array<T, 3> gradient((T)0, (T)0, (T)0);
        for (plint iPop = 1; iPop < D::q; ++iPop) {
            plint nextX = iX + D::c[iPop][0];
            plint nextY = iY + D::c[iPop][1];
            plint nextZ = iZ + D::c[iPop][2];

            T nextVal = 0.0;

            // First, extrapolate the scalar field on the wall (if necessary).
            if (isAnyWall(flag_->get(
                    nextX + relativeOffsetFS.x, nextY + relativeOffsetFS.y,
                    nextZ + relativeOffsetFS.z)))
            {
                T locVal = scalar.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z);
                plint opp = indexTemplates::opposite<D>(iPop);
                plint prevX = iX + D::c[opp][0];
                plint prevY = iY + D::c[opp][1];
                plint prevZ = iZ + D::c[opp][2];
                if (isAnyWall(flag_->get(
                        prevX + relativeOffsetFS.x, prevY + relativeOffsetFS.y,
                        prevZ + relativeOffsetFS.z)))
                {
                    nextVal = locVal;
                } else {
                    T prevVal = scalar.get(prevX + ofs.x, prevY + ofs.y, prevZ + ofs.z);
                    nextVal = (T)2 * locVal - prevVal;
                }
            } else {
                nextVal = scalar.get(nextX + ofs.x, nextY + ofs.y, nextZ + ofs.z);
            }

            gradient[0] += D::t[iPop] * D::c[iPop][0] * nextVal;
            gradient[1] += D::t[iPop] * D::c[iPop][1] * nextVal;
            gradient[2] += D::t[iPop] * D::c[iPop][2] * nextVal;
        }
        gradient *= D::invCs2;

        return (gradient);
    }

    // Compute the divergence of a tensor field "the lattice Boltzmann way".
    // With this method the tensor field values are first extrapolated on the neighboring wall cells
    // (if necessary). So, in contrast to the above finite difference implementation, the wall cells
    // are not excluded here. The following function can be called with a lattice descriptor
    // "Descriptor2", which is different from the lattice descriptor "Descriptor" used in the
    // simulation.
    template <typename T2, template <typename U2> class Descriptor2>
    T lbmComputeDivergence(
        TensorField3D<T, 3> const &tensor, Dot3D const &ofs, plint iX, plint iY, plint iZ)
    {
        typedef Descriptor2<T2> D;
        using namespace freeSurfaceFlag;

        // PLB_ASSERT(!isAnyWall(flag_->get(iX + relativeOffsetFS.x, iY + relativeOffsetFS.y, iZ +
        // relativeOffsetFS.z)));

        T divergence = 0.0;
        for (plint iPop = 1; iPop < D::q; ++iPop) {
            plint nextX = iX + D::c[iPop][0];
            plint nextY = iY + D::c[iPop][1];
            plint nextZ = iZ + D::c[iPop][2];

            Array<T, 3> nextVal((T)0, (T)0, (T)0);

            // First, extrapolate the scalar field on the wall (if necessary).
            if (isAnyWall(flag_->get(
                    nextX + relativeOffsetFS.x, nextY + relativeOffsetFS.y,
                    nextZ + relativeOffsetFS.z)))
            {
                Array<T, 3> locVal = tensor.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z);
                plint opp = indexTemplates::opposite<D>(iPop);
                plint prevX = iX + D::c[opp][0];
                plint prevY = iY + D::c[opp][1];
                plint prevZ = iZ + D::c[opp][2];
                if (isAnyWall(flag_->get(
                        prevX + relativeOffsetFS.x, prevY + relativeOffsetFS.y,
                        prevZ + relativeOffsetFS.z)))
                {
                    nextVal = locVal;
                } else {
                    Array<T, 3> prevVal = tensor.get(prevX + ofs.x, prevY + ofs.y, prevZ + ofs.z);
                    nextVal = (T)2 * locVal - prevVal;
                }
            } else {
                nextVal = tensor.get(nextX + ofs.x, nextY + ofs.y, nextZ + ofs.z);
            }

            divergence += D::t[iPop]
                          * (D::c[iPop][0] * nextVal[0] + D::c[iPop][1] * nextVal[1]
                             + D::c[iPop][2] * nextVal[2]);
        }
        divergence *= D::invCs2;

        return (divergence);
    }

    std::map<Node, T> &massExcess()
    {
        return interfaceLists_->massExcess;
    }
    std::set<Node> &interfaceToFluid()
    {
        return interfaceLists_->interfaceToFluid;
    }
    std::set<Node> &interfaceToEmpty()
    {
        return interfaceLists_->interfaceToEmpty;
    }
    std::set<Node> &emptyToInterface()
    {
        return interfaceLists_->emptyToInterface;
    }

    BlockLattice3D<T, Descriptor> *latticeP()
    {
        return fluid_;
    }
    ScalarField3D<T> *rhoBarP()
    {
        return rhoBar_;
    }
    TensorField3D<T, 3> *jP()
    {
        return j_;
    }
    ScalarField3D<T> *massP()
    {
        return mass_;
    }
    ScalarField3D<T> *volumeFractionP()
    {
        return volumeFraction_;
    }
    ScalarField3D<int> *flagP()
    {
        return flag_;
    }
    TensorField3D<T, 3> *normalP()
    {
        return normal_;
    }
    ScalarField3D<T> *curvatureP()
    {
        return curvature_;
    }
    ScalarField3D<T> *outsideDensityP()
    {
        return outsideDensity_;
    }

    Dot3D const &absOffset() const
    {
        return absoluteOffset;
    }
    Dot3D const &rhoBarOffset() const
    {
        return relativeOffsetRhoBar;
    }
    Dot3D const &jOffset() const
    {
        return relativeOffsetJ;
    }
    Dot3D const &massOffset() const
    {
        return relativeOffsetMass;
    }
    Dot3D const &volumeFractionOffset() const
    {
        return relativeOffsetVF;
    }
    Dot3D const &flagOffset() const
    {
        return relativeOffsetFS;
    }
    Dot3D const &normalOffset() const
    {
        return relativeOffsetNormal;
    }
    Dot3D const &curvatureOffset() const
    {
        return relativeOffsetC;
    }
    Dot3D const &outsideDensityOffset() const
    {
        return relativeOffsetOD;
    }

    Box3D getBoundingBox() const
    {
        return volumeFraction_->getBoundingBox();
    }

private:
    BlockLattice3D<T, Descriptor> *fluid_;
    ScalarField3D<T> *rhoBar_;
    TensorField3D<T, 3> *j_;
    ScalarField3D<T> *mass_;
    ScalarField3D<T> *volumeFraction_;
    ScalarField3D<int> *flag_;
    TensorField3D<T, 3> *normal_;
    AtomicContainerBlock3D *containerInterfaceLists_;
    InterfaceLists<T, Descriptor> *interfaceLists_;
    ScalarField3D<T> *curvature_;
    ScalarField3D<T> *outsideDensity_;

    Dot3D absoluteOffset, relativeOffsetRhoBar, relativeOffsetJ, relativeOffsetMass,
        relativeOffsetVF, relativeOffsetFS, relativeOffsetNormal, relativeOffsetC, relativeOffsetOD;
};

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeNormalizedVolumeFraction(
    MultiScalarField3D<T> &volumeFraction, Box3D domain)
{
    std::unique_ptr<MultiScalarField3D<T> > result = extractSubDomain<T>(volumeFraction, domain);
    boundScalarField(*result, domain, (T)0, (T)1);
    return result;
}

template <typename T>
std::unique_ptr<MultiScalarField3D<T> > computeNormalizedVolumeFraction(
    MultiScalarField3D<T> &volumeFraction)
{
    return computeNormalizedVolumeFraction(volumeFraction, volumeFraction.getBoundingBox());
}

}  // namespace plb

#endif  // FREE_SURFACE_UTIL_3D_H
