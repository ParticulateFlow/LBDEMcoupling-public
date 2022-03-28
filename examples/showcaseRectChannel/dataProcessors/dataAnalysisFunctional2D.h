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
#ifndef DATA_ANALYSIS_FUNCTIONAL_2D_H
#define DATA_ANALYSIS_FUNCTIONAL_2D_H

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for BlockLattice ******* */

template <typename T, template <typename U> class Descriptor>
class BoxSumRhoBarFunctional2D : public ReductiveBoxProcessingFunctional2D_L<T, Descriptor> {
public:
    BoxSumRhoBarFunctional2D();
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual BoxSumRhoBarFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    T getSumRhoBar() const;

private:
    plint sumRhoBarId;
};

template <typename T, template <typename U> class Descriptor>
class DotSumRhoBarFunctional2D : public ReductiveDotProcessingFunctional2D_L<T, Descriptor> {
public:
    DotSumRhoBarFunctional2D();
    virtual void process(DotList2D const &dotList, BlockLattice2D<T, Descriptor> &lattice);
    virtual DotSumRhoBarFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    T getSumRhoBar() const;

private:
    plint sumRhoBarId;
};

template <typename T, template <typename U> class Descriptor>
class BoxSumEnergyFunctional2D : public ReductiveBoxProcessingFunctional2D_L<T, Descriptor> {
public:
    BoxSumEnergyFunctional2D();
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual BoxSumEnergyFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    T getSumEnergy() const;
    virtual void getDimensionsX(std::vector<int> &dimensions) const;
    virtual void getDimensionsT(std::vector<int> &dimensions) const;

private:
    plint sumEnergyId;
};

template <typename T, template <typename U> class Descriptor>
class DotSumVelocityComponentFunctional2D :
    public ReductiveDotProcessingFunctional2D_L<T, Descriptor> {
public:
    DotSumVelocityComponentFunctional2D(int iComponent_);
    virtual void process(DotList2D const &dotList, BlockLattice2D<T, Descriptor> &lattice);
    virtual DotSumVelocityComponentFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    T getSumVelocityComponent() const;

private:
    plint sumVelocityComponentId;
    int iComponent;
};

template <typename T, template <typename U> class Descriptor, class BoolMask>
class CountLatticeElementsFunctional2D :
    public ReductiveBoxProcessingFunctional2D_L<T, Descriptor> {
public:
    CountLatticeElementsFunctional2D(BoolMask boolMask_);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    virtual CountLatticeElementsFunctional2D<T, Descriptor, BoolMask> *clone() const;
    plint getCount() const;

private:
    plint countId;
    BoolMask boolMask;
};

/* *************** Data Functionals for BlockLattice ***************** */

template <typename T, template <typename U> class Descriptor>
class CopyPopulationsFunctional2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor, T, Descriptor> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &latticeFrom,
        BlockLattice2D<T, Descriptor> &latticeTo);
    virtual CopyPopulationsFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <
    typename T1, template <typename U1> class Descriptor1, typename T2,
    template <typename U2> class Descriptor2>
class CopyConvertPopulationsFunctional2D :
    public BoxProcessingFunctional2D_LL<T1, Descriptor1, T2, Descriptor2> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T1, Descriptor1> &latticeFrom,
        BlockLattice2D<T2, Descriptor2> &latticeTo);
    virtual CopyConvertPopulationsFunctional2D<T1, Descriptor1, T2, Descriptor2> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, template <typename U> class Descriptor>
class LatticeCopyAllFunctional2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor, T, Descriptor> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &latticeFrom,
        BlockLattice2D<T, Descriptor> &latticeTo);
    virtual LatticeCopyAllFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class LatticeRegenerateFunctional2D :
    public BoxProcessingFunctional2D_LL<T, Descriptor, T, Descriptor> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &latticeFrom,
        BlockLattice2D<T, Descriptor> &latticeTo);
    virtual LatticeRegenerateFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, template <typename U> class Descriptor>
class BoxDensityFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField);
    virtual BoxDensityFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxRhoBarFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField);
    virtual BoxRhoBarFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxRhoBarJfunctional2D : public BoxProcessingFunctional2D {
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> fields);
    virtual BoxRhoBarJfunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, template <typename U> class Descriptor>
class PackedRhoBarJfunctional2D : public BoxProcessingFunctional2D_LN<T, Descriptor, T> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, NTensorField2D<T> &rhoBarJ);
    virtual PackedRhoBarJfunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, template <typename U> class Descriptor>
class BoxKineticEnergyFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField);
    virtual BoxKineticEnergyFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxVelocityNormFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField);
    virtual BoxVelocityNormFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxVelocityComponentFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    BoxVelocityComponentFunctional2D(int iComponent_);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField);
    virtual BoxVelocityComponentFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    int iComponent;
};

template <typename T, template <typename U> class Descriptor>
class BoxVelocityFunctional2D :
    public BoxProcessingFunctional2D_LT<T, Descriptor, T, Descriptor<T>::d> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        TensorField2D<T, Descriptor<T>::d> &tensorField);
    virtual BoxVelocityFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxPiNeqFunctional2D :
    public BoxProcessingFunctional2D_LT<T, Descriptor, T, SymmetricTensor<T, Descriptor>::n> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq);
    virtual BoxPiNeqFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxShearStressFunctional2D :
    public BoxProcessingFunctional2D_LT<T, Descriptor, T, SymmetricTensor<T, Descriptor>::n> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &PiNeq);
    virtual BoxShearStressFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxStrainRateFromStressFunctional2D :
    public BoxProcessingFunctional2D_LT<T, Descriptor, T, SymmetricTensor<T, Descriptor>::n> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        TensorField2D<T, SymmetricTensor<T, Descriptor>::n> &S);
    virtual BoxStrainRateFromStressFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxTemperatureFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField);
    virtual BoxTemperatureFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxSoundSpeedFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField);
    virtual BoxSoundSpeedFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxPopulationFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    BoxPopulationFunctional2D(plint iComponent_);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &population);
    virtual BoxPopulationFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    plint iComponent;
};

template <typename T, template <typename U> class Descriptor>
class BoxEquilibriumFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    BoxEquilibriumFunctional2D(plint iComponent_);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &equilibrium);
    virtual BoxEquilibriumFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    plint iComponent;
};

template <typename T, template <typename U> class Descriptor>
class BoxAllPopulationsFunctional2D :
    public BoxProcessingFunctional2D_LT<T, Descriptor, T, Descriptor<T>::q> {
public:
    BoxAllPopulationsFunctional2D();
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        TensorField2D<T, Descriptor<T>::q> &populations);
    virtual BoxAllPopulationsFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxAllPopulationsToLatticeFunctional2D :
    public BoxProcessingFunctional2D_LT<T, Descriptor, T, Descriptor<T>::q> {
public:
    BoxAllPopulationsToLatticeFunctional2D();
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        TensorField2D<T, Descriptor<T>::q> &populations);
    virtual BoxAllPopulationsToLatticeFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxOmegaFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField);
    virtual BoxOmegaFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, template <typename U> class Descriptor>
class BoxExternalForceFunctional2D :
    public BoxProcessingFunctional2D_LT<T, Descriptor, T, Descriptor<T>::d> {
public:
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        TensorField2D<T, Descriptor<T>::d> &tensorField);
    virtual BoxExternalForceFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T, template <typename U> class Descriptor>
class BoxExternalScalarFunctional2D : public BoxProcessingFunctional2D_LS<T, Descriptor, T> {
public:
    BoxExternalScalarFunctional2D(int whichScalar_);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice, ScalarField2D<T> &scalarField);
    virtual BoxExternalScalarFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    int whichScalar;
};

template <typename T, template <typename U> class Descriptor>
class BoxExternalVectorFunctional2D :
    public BoxProcessingFunctional2D_LT<T, Descriptor, T, Descriptor<T>::d> {
public:
    BoxExternalVectorFunctional2D(int vectorBeginsAt_);
    virtual void process(
        Box2D domain, BlockLattice2D<T, Descriptor> &lattice,
        TensorField2D<T, Descriptor<T>::d> &tensorField);
    virtual BoxExternalVectorFunctional2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    int vectorBeginsAt;
};

/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for ScalarField ******** */

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the sum is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template <typename T>
class BoxScalarSumFunctional2D : public ReductiveBoxProcessingFunctional2D_S<T> {
public:
    BoxScalarSumFunctional2D();
    virtual void process(Box2D domain, ScalarField2D<T> &scalarField);
    virtual BoxScalarSumFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    T getSumScalar() const;

private:
    plint sumScalarId;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the average is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template <typename T>
class MaskedBoxScalarAverageFunctional2D : public ReductiveBoxProcessingFunctional2D_SS<T, int> {
public:
    MaskedBoxScalarAverageFunctional2D(int flag_);
    virtual void process(Box2D domain, ScalarField2D<T> &scalarField, ScalarField2D<int> &mask);
    virtual MaskedBoxScalarAverageFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getAverageScalar() const;

private:
    plint averageScalarId;
    int flag;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the min is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template <typename T>
class BoxScalarMinFunctional2D : public ReductiveBoxProcessingFunctional2D_S<T> {
public:
    BoxScalarMinFunctional2D();
    virtual void process(Box2D domain, ScalarField2D<T> &scalarField);
    virtual BoxScalarMinFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    T getMinScalar() const;

private:
    plint maxScalarId;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the max is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template <typename T>
class MaskedBoxScalarMinFunctional2D : public ReductiveBoxProcessingFunctional2D_SS<T, int> {
public:
    MaskedBoxScalarMinFunctional2D(int flag_);
    virtual void process(Box2D domain, ScalarField2D<T> &scalarField, ScalarField2D<int> &mask);
    virtual MaskedBoxScalarMinFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getMinScalar() const;

private:
    plint maxScalarId;
    int flag;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the max is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template <typename T>
class BoxScalarMaxFunctional2D : public ReductiveBoxProcessingFunctional2D_S<T> {
public:
    BoxScalarMaxFunctional2D();
    virtual void process(Box2D domain, ScalarField2D<T> &scalarField);
    virtual BoxScalarMaxFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    T getMaxScalar() const;

private:
    plint maxScalarId;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the max is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template <typename T>
class MaskedBoxScalarMaxFunctional2D : public ReductiveBoxProcessingFunctional2D_SS<T, int> {
public:
    MaskedBoxScalarMaxFunctional2D(int flag_);
    virtual void process(Box2D domain, ScalarField2D<T> &scalarField, ScalarField2D<int> &mask);
    virtual MaskedBoxScalarMaxFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getMaxScalar() const;

private:
    plint maxScalarId;
    int flag;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the sum is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template <typename T>
class BoundedBoxScalarSumFunctional2D : public BoundedReductiveBoxProcessingFunctional2D_S<T> {
public:
    BoundedBoxScalarSumFunctional2D();
    virtual void processBulk(Box2D domain, ScalarField2D<T> &scalarField);
    virtual void processEdge(
        int direction, int orientation, Box2D domain, ScalarField2D<T> &scalarField);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, ScalarField2D<T> &scalarField);
    virtual BoundedBoxScalarSumFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    T getSumScalar() const;

private:
    plint sumScalarId;
};

template <typename T, class BoolMask>
class CountScalarElementsFunctional2D : public ReductiveBoxProcessingFunctional2D_S<T> {
public:
    CountScalarElementsFunctional2D(BoolMask boolMask_);
    virtual void process(Box2D domain, ScalarField2D<T> &field);
    virtual CountScalarElementsFunctional2D<T, BoolMask> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    plint getCount() const;

private:
    plint countId;
    BoolMask boolMask;
};

/* *************** Data Functionals for scalar-fields **************** */

template <typename T>
class ComputeAbsoluteValueFunctional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B);
    virtual ComputeAbsoluteValueFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class ComputeScalarSqrtFunctional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B);
    virtual ComputeScalarSqrtFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class ComputeScalarPowFunctional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    ComputeScalarPowFunctional2D(T power_);
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B);
    virtual ComputeScalarPowFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T power;
};

template <typename T>
class ComputeScalarLogFunctional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B);
    virtual ComputeScalarLogFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T1, typename T2>
class CopyConvertScalarFunctional2D : public BoxProcessingFunctional2D_SS<T1, T2> {
public:
    virtual void process(Box2D domain, ScalarField2D<T1> &field1, ScalarField2D<T2> &field2);
    virtual CopyConvertScalarFunctional2D<T1, T2> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class ExtractScalarSubDomainFunctional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void process(Box2D domain, ScalarField2D<T> &field1, ScalarField2D<T> &field2);
    virtual ExtractScalarSubDomainFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, class Function>
class ApplyScalarFunctional2D : public BoxProcessingFunctional2D_S<T> {
public:
    ApplyScalarFunctional2D(Function f_);
    virtual void process(Box2D domain, ScalarField2D<T> &field);
    virtual ApplyScalarFunctional2D<T, Function> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    Function f;
};

template <typename T, class EvalFunction>
class EvaluateScalarFunctional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    EvaluateScalarFunctional2D(EvalFunction f_);
    virtual void process(Box2D domain, ScalarField2D<T> &field, ScalarField2D<T> &result);
    virtual EvaluateScalarFunctional2D<T, EvalFunction> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    EvalFunction f;
};

template <typename T>
class A_lt_alpha_functional2D : public BoxProcessingFunctional2D_SS<T, int> {
public:
    A_lt_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<int> &result);
    virtual A_lt_alpha_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class A_gt_alpha_functional2D : public BoxProcessingFunctional2D_SS<T, int> {
public:
    A_gt_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<int> &result);
    virtual A_gt_alpha_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class A_plus_alpha_functional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    A_plus_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result);
    virtual A_plus_alpha_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class A_minus_alpha_functional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    A_minus_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result);
    virtual A_minus_alpha_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class Alpha_minus_A_functional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    Alpha_minus_A_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result);
    virtual Alpha_minus_A_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class A_times_alpha_functional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    A_times_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result);
    virtual A_times_alpha_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class A_dividedBy_alpha_functional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    A_dividedBy_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result);
    virtual A_dividedBy_alpha_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class Alpha_dividedBy_A_functional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    Alpha_dividedBy_A_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &result);
    virtual Alpha_dividedBy_A_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class A_plus_alpha_inplace_functional2D : public BoxProcessingFunctional2D_S<T> {
public:
    A_plus_alpha_inplace_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A);
    virtual A_plus_alpha_inplace_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class A_minus_alpha_inplace_functional2D : public BoxProcessingFunctional2D_S<T> {
public:
    A_minus_alpha_inplace_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A);
    virtual A_minus_alpha_inplace_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class A_times_alpha_inplace_functional2D : public BoxProcessingFunctional2D_S<T> {
public:
    A_times_alpha_inplace_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A);
    virtual A_times_alpha_inplace_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class A_dividedBy_alpha_inplace_functional2D : public BoxProcessingFunctional2D_S<T> {
public:
    A_dividedBy_alpha_inplace_functional2D(T alpha_);
    virtual void process(Box2D domain, ScalarField2D<T> &A);
    virtual A_dividedBy_alpha_inplace_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T>
class A_lt_B_functional2D : public BoxProcessingFunctional2D {
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> scalarFields);
    virtual A_lt_B_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class A_gt_B_functional2D : public BoxProcessingFunctional2D {
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D *> scalarFields);
    virtual A_gt_B_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class A_plus_B_functional2D : public ScalarFieldBoxProcessingFunctional2D<T> {
public:
    virtual void process(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual A_plus_B_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class A_minus_B_functional2D : public ScalarFieldBoxProcessingFunctional2D<T> {
public:
    virtual void process(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual A_minus_B_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class A_times_B_functional2D : public ScalarFieldBoxProcessingFunctional2D<T> {
public:
    virtual void process(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual A_times_B_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class A_dividedBy_B_functional2D : public ScalarFieldBoxProcessingFunctional2D<T> {
public:
    virtual void process(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual A_dividedBy_B_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class A_plus_B_inplace_functional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B);
    virtual A_plus_B_inplace_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class A_minus_B_inplace_functional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B);
    virtual A_minus_B_inplace_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class A_times_B_inplace_functional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B);
    virtual A_times_B_inplace_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class A_dividedBy_B_inplace_functional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void process(Box2D domain, ScalarField2D<T> &A, ScalarField2D<T> &B);
    virtual A_dividedBy_B_inplace_functional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

template <typename T, int nDim>
class ComputeTensorSqrtFunctional2D : public BoxProcessingFunctional2D_TT<T, nDim, T, nDim> {
public:
    virtual void process(Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);
    virtual ComputeTensorSqrtFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim, class BoolMask>
class CountTensorElementsFunctional2D : public ReductiveBoxProcessingFunctional2D_T<T, nDim> {
public:
    CountTensorElementsFunctional2D(BoolMask boolMask_);
    virtual void process(Box2D domain, TensorField2D<T, nDim> &field);
    virtual CountTensorElementsFunctional2D<T, nDim, BoolMask> *clone() const;
    plint getCount() const;

private:
    plint countId;
    BoolMask boolMask;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the sum is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template <typename T, int nDim>
class BoxTensorSumFunctional2D : public ReductiveBoxProcessingFunctional2D_T<T, nDim> {
public:
    BoxTensorSumFunctional2D();
    virtual void process(Box2D domain, TensorField2D<T, nDim> &tensorField);
    virtual BoxTensorSumFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
    }
    Array<T, nDim> getSumTensor() const;

private:
    Array<plint, nDim> sumTensorId;
};

template <typename T1, typename T2, int nDim>
class CopyConvertTensorFunctional2D : public BoxProcessingFunctional2D_TT<T1, nDim, T2, nDim> {
public:
    virtual void process(
        Box2D domain, TensorField2D<T1, nDim> &field1, TensorField2D<T2, nDim> &field2);
    virtual CopyConvertTensorFunctional2D<T1, T2, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class ExtractTensorSubDomainFunctional2D : public BoxProcessingFunctional2D_TT<T, nDim, T, nDim> {
public:
    virtual void process(
        Box2D domain, TensorField2D<T, nDim> &field1, TensorField2D<T, nDim> &field2);
    virtual ExtractTensorSubDomainFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class ExtractTensorComponentFunctional2D : public BoxProcessingFunctional2D_ST<T, T, nDim> {
public:
    ExtractTensorComponentFunctional2D(int iComponent_);
    virtual void process(
        Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, nDim> &tensorField);
    virtual ExtractTensorComponentFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    int iComponent;
};

template <typename T, int nDim>
class ComputeNormFunctional2D : public BoxProcessingFunctional2D_ST<T, T, nDim> {
public:
    virtual void process(
        Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, nDim> &tensorField);
    virtual ComputeNormFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class ComputeNormSqrFunctional2D : public BoxProcessingFunctional2D_ST<T, T, nDim> {
public:
    virtual void process(
        Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, nDim> &tensorField);
    virtual ComputeNormSqrFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class ComputeSymmetricTensorNormFunctional2D : public BoxProcessingFunctional2D_ST<T, T, 3> {
public:
    virtual void process(
        Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, 3> &tensorField);
    virtual ComputeSymmetricTensorNormFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class ComputeSymmetricTensorNormSqrFunctional2D : public BoxProcessingFunctional2D_ST<T, T, 3> {
public:
    virtual void process(
        Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, 3> &tensorField);
    virtual ComputeSymmetricTensorNormSqrFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class ComputeSymmetricTensorTraceFunctional2D : public BoxProcessingFunctional2D_ST<T, T, 3> {
public:
    virtual void process(
        Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, 3> &tensorField);
    virtual ComputeSymmetricTensorTraceFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class BoxLocalMaximumPerComponentFunctional2D : public BoxProcessingFunctional2D_ST<T, T, nDim> {
public:
    virtual void process(
        Box2D domain, ScalarField2D<T> &scalarField, TensorField2D<T, nDim> &tensorField);
    virtual BoxLocalMaximumPerComponentFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T>
class BoxGradientFunctional2D : public BoundedBoxProcessingFunctional2D_ST<T, T, 2> {
public:
    virtual void processBulk(Box2D domain, ScalarField2D<T> &phi, TensorField2D<T, 2> &gradient);

    virtual void processEdge(
        int direction, int orientation, Box2D domain, ScalarField2D<T> &phi,
        TensorField2D<T, 2> &gradient);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, ScalarField2D<T> &phi,
        TensorField2D<T, 2> &gradient);

    virtual BoxGradientFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class BoxBulkVorticityFunctional2D : public BoxProcessingFunctional2D_ST<T, T, nDim> {
public:
    virtual void process(
        Box2D domain, ScalarField2D<T> &vorticity, TensorField2D<T, nDim> &velocity);
    virtual BoxBulkVorticityFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class BoxVorticityFunctional2D : public BoundedBoxProcessingFunctional2D_ST<T, T, nDim> {
public:
    virtual void processBulk(
        Box2D domain, ScalarField2D<T> &vorticity, TensorField2D<T, nDim> &velocity);
    virtual void processEdge(
        int direction, int orientation, Box2D domain, ScalarField2D<T> &vorticity,
        TensorField2D<T, nDim> &velocity);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, ScalarField2D<T> &vorticity,
        TensorField2D<T, nDim> &velocity);
    virtual BoxVorticityFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class BoxBulkStrainRateFunctional2D :
    public BoxProcessingFunctional2D_TT<T, nDim, T, SymmetricTensorImpl<T, nDim>::n> {
public:
    virtual void process(
        Box2D domain, TensorField2D<T, nDim> &velocity,
        TensorField2D<T, SymmetricTensorImpl<T, nDim>::n> &S);
    virtual BoxBulkStrainRateFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class BoxStrainRateFunctional2D :
    public BoundedBoxProcessingFunctional2D_TT<T, nDim, T, SymmetricTensorImpl<T, nDim>::n> {
public:
    virtual void processBulk(
        Box2D domain, TensorField2D<T, nDim> &velocity,
        TensorField2D<T, SymmetricTensorImpl<T, nDim>::n> &S);
    virtual void processEdge(
        int direction, int orientation, Box2D domain, TensorField2D<T, nDim> &velocity,
        TensorField2D<T, SymmetricTensorImpl<T, nDim>::n> &S);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, TensorField2D<T, nDim> &velocity,
        TensorField2D<T, SymmetricTensorImpl<T, nDim>::n> &S);
    virtual BoxStrainRateFunctional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class Tensor_A_times_alpha_functional2D : public BoxProcessingFunctional2D_TT<T, nDim, T, nDim> {
public:
    Tensor_A_times_alpha_functional2D(T alpha_);
    virtual void process(Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &result);
    virtual Tensor_A_times_alpha_functional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T, int nDim>
class Tensor_A_plus_B_functional2D : public TensorFieldBoxProcessingFunctional2D<T, nDim> {
public:
    virtual void process(Box2D domain, std::vector<TensorField2D<T, nDim> *> tensorFields);
    virtual Tensor_A_plus_B_functional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class Tensor_A_minus_B_functional2D : public TensorFieldBoxProcessingFunctional2D<T, nDim> {
public:
    virtual void process(Box2D domain, std::vector<TensorField2D<T, nDim> *> tensorFields);
    virtual Tensor_A_minus_B_functional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class Tensor_A_times_B_functional2D : public TensorFieldBoxProcessingFunctional2D<T, nDim> {
public:
    virtual void process(Box2D domain, std::vector<TensorField2D<T, nDim> *> tensorFields);
    virtual Tensor_A_times_B_functional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class Tensor_A_dividedBy_B_functional2D : public TensorFieldBoxProcessingFunctional2D<T, nDim> {
public:
    virtual void process(Box2D domain, std::vector<TensorField2D<T, nDim> *> tensorFields);
    virtual Tensor_A_dividedBy_B_functional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class Tensor_A_plus_B_inplace_functional2D : public BoxProcessingFunctional2D_TT<T, nDim, T, nDim> {
public:
    virtual void process(Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);
    virtual Tensor_A_plus_B_inplace_functional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class Tensor_A_minus_B_inplace_functional2D :
    public BoxProcessingFunctional2D_TT<T, nDim, T, nDim> {
public:
    virtual void process(Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);
    virtual Tensor_A_minus_B_inplace_functional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class Tensor_A_times_B_inplace_functional2D :
    public BoxProcessingFunctional2D_TT<T, nDim, T, nDim> {
public:
    virtual void process(Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);
    virtual Tensor_A_times_B_inplace_functional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T, int nDim>
class Tensor_A_times_alpha_inplace_functional2D : public BoxProcessingFunctional2D_T<T, nDim> {
public:
    Tensor_A_times_alpha_inplace_functional2D(T alpha_);
    virtual void process(Box2D domain, TensorField2D<T, nDim> &A);
    virtual Tensor_A_times_alpha_inplace_functional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T alpha;
};

template <typename T, int nDim>
class Tensor_A_dividedBy_B_inplace_functional2D :
    public BoxProcessingFunctional2D_TT<T, nDim, T, nDim> {
public:
    virtual void process(Box2D domain, TensorField2D<T, nDim> &A, TensorField2D<T, nDim> &B);
    virtual Tensor_A_dividedBy_B_inplace_functional2D<T, nDim> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONAL_2D_H

// Explicitly include generic algorithms which are never precompiled (not even in precompiled
// version)
#include "dataProcessors/dataAnalysisGenerics2D.h"
