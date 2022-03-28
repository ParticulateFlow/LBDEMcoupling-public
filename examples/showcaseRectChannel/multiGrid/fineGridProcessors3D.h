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
 * Dynamics and data processors used to implement 3D grid refinement -- header file.
 */

#ifndef FINE_GRID_PROCESSORS_3D_H
#define FINE_GRID_PROCESSORS_3D_H

#include <vector>

#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiGrid/gridRefinement.h"
#include "multiGrid/gridRefinementDynamics.h"

namespace plb {

/// Processor to be added to coarse lattice: copies t1 values to t0 after each numTimeSteps
/// iteration
template <typename T, template <typename U> class Descriptor>
class Copy_t1_to_t0_3D : public BoxProcessingFunctional3D_L<T, Descriptor> {
public:
    Copy_t1_to_t0_3D(plint numTimeSteps_, plint executionTime_);
    virtual void process(Box3D domain, BlockLattice3D<T, Descriptor> &fineLattice);
    virtual Copy_t1_to_t0_3D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::dataStructure;
    }

private:
    plint numTimeSteps;
    plint executionTime;
};

//////////////////////////////////////////////////////////////////////////////
/////////////////////          YZ plane           ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/// Bulk interpolation
template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationYZ :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationYZ(RescaleEngine<T, Descriptor> *rescaleEngine_) :
        rescaleEngine(rescaleEngine_)
    { }

    virtual ~ScalarCubicInterpolationYZ()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationYZ<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

/// The edges
template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationYZLineY3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationYZLineY3D(plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_);

    virtual ~ScalarCubicInterpolationYZLineY3D()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationYZLineY3D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
    plint delta;
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationYZLineZ3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationYZLineZ3D(plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_);

    virtual ~ScalarCubicInterpolationYZLineZ3D()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationYZLineZ3D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
private:
    plint delta;
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

/// The corners
template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationYZCorner3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationYZCorner3D(
        plint deltaY_, plint deltaZ_, RescaleEngine<T, Descriptor> *rescaleEngine_);

    virtual ~ScalarCubicInterpolationYZCorner3D()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationYZCorner3D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
    plint deltaY;
    plint deltaZ;
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

//////////////////////////////////////////////////////////////////////////////
/////////////////////          XY plane           ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/// Bulk interpolation
template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationXY :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationXY(RescaleEngine<T, Descriptor> *rescaleEngine_) :
        rescaleEngine(rescaleEngine_)
    { }

    virtual ~ScalarCubicInterpolationXY()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationXY<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

/// The edges
template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationXYLineX3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationXYLineX3D(plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_);

    virtual ~ScalarCubicInterpolationXYLineX3D()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationXYLineX3D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
    plint delta;
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationXYLineY3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationXYLineY3D(plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_);

    virtual ~ScalarCubicInterpolationXYLineY3D()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationXYLineY3D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
private:
    plint delta;
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

/// The corners
template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationXYCorner3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationXYCorner3D(
        plint deltaX_, plint deltaY_, RescaleEngine<T, Descriptor> *rescaleEngine_);

    virtual ~ScalarCubicInterpolationXYCorner3D()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationXYCorner3D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
    plint deltaX;
    plint deltaY;
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

//////////////////////////////////////////////////////////////////////////////
/////////////////////          XZ plane           ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/// Bulk interpolation
template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationXZ :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationXZ(RescaleEngine<T, Descriptor> *rescaleEngine_) :
        rescaleEngine(rescaleEngine_)
    { }

    virtual ~ScalarCubicInterpolationXZ()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationXZ<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

/// The edges
template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationXZLineX3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationXZLineX3D(plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_);

    virtual ~ScalarCubicInterpolationXZLineX3D()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationXZLineX3D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
    plint delta;
    RescaleEngine<T, Descriptor> *rescaleEngine;

    Array<Array<T, 2>, 4> points;  // where the interpolated values will be located
};

template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationXZLineZ3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationXZLineZ3D(plint delta_, RescaleEngine<T, Descriptor> *rescaleEngine_);

    virtual ~ScalarCubicInterpolationXZLineZ3D()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationXZLineZ3D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
private:
    plint delta;
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

/// The corners
template <typename T, template <typename U> class Descriptor>
class ScalarCubicInterpolationXZCorner3D :
    public BoxProcessingFunctional3D_LL<T, Descriptor, T, Descriptor> {
public:
    ScalarCubicInterpolationXZCorner3D(
        plint deltaX_, plint deltaZ_, RescaleEngine<T, Descriptor> *rescaleEngine_);

    virtual ~ScalarCubicInterpolationXZCorner3D()
    {
        delete rescaleEngine;
    }

    virtual void process(
        Box3D coarseDomain, BlockLattice3D<T, Descriptor> &coarseLattice,
        BlockLattice3D<T, Descriptor> &fineLattice);

    virtual ScalarCubicInterpolationXZCorner3D<T, Descriptor> *clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::dataStructure;
    }

private:
    plint deltaX;
    plint deltaZ;
    RescaleEngine<T, Descriptor> *rescaleEngine;
};

}  // namespace plb

#endif  // FINE_GRID_PROCESSORS_3D_H
