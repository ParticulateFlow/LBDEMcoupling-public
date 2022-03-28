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

#ifndef DYNAMICS_PROCESSOR_3D_H
#define DYNAMICS_PROCESSOR_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

/* *************** Class ExternalRhoJcollideAndStream3D ******************* */

template <typename T, template <typename U> class Descriptor>
class ExternalRhoJcollideAndStream3D : public BoxProcessingFunctional3D {
public:
    // Block 0: lattice; Block 1: rhoBar; Block 2: j.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual ExternalRhoJcollideAndStream3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    void collide(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
        ScalarField3D<T> const &rhoBarField, Dot3D const &offset1,
        TensorField3D<T, 3> const &jField, Dot3D const &offset2, BlockStatistics &stat);
    void bulkCollideAndStream(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
        ScalarField3D<T> const &rhoBarField, Dot3D const &offset1,
        TensorField3D<T, 3> const &jField, Dot3D const &offset2, BlockStatistics &stat);
    void boundaryStream(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &bound, Box3D const &domain);
};

/* *************** Class PackedExternalRhoJcollideAndStream3D ******************* */

template <typename T, template <typename U> class Descriptor>
class PackedExternalRhoJcollideAndStream3D : public BoxProcessingFunctional3D_LN<T, Descriptor, T> {
public:
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, NTensorField3D<T> &rhoBarJfield);
    virtual PackedExternalRhoJcollideAndStream3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    void collide(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
        NTensorField3D<T> const &rhoBarJfield, Dot3D const &offset, BlockStatistics &stat);
    void bulkCollideAndStream(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
        NTensorField3D<T> const &rhoBarJfield, Dot3D const &offset, BlockStatistics &stat);
    void boundaryStream(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &bound, Box3D const &domain);
};

/* *************** Class WaveAbsorptionExternalRhoJcollideAndStream3D ******************* */

template <typename T, template <typename U> class Descriptor>
class WaveAbsorptionExternalRhoJcollideAndStream3D : public BoxProcessingFunctional3D {
public:
    WaveAbsorptionExternalRhoJcollideAndStream3D(T rhoBarF_, Array<T, Descriptor<T>::d> const &uF_);
    // Block 0: lattice; Block 1: rhoBar; Block 2: j; Block 3: sigma.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual WaveAbsorptionExternalRhoJcollideAndStream3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    void collide(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
        ScalarField3D<T> const &rhoBarField, Dot3D const &offset1,
        TensorField3D<T, Descriptor<T>::d> const &jField, Dot3D const &offset2,
        ScalarField3D<T> const &sigmaField, Dot3D const &offset3, BlockStatistics &stat);
    void bulkCollideAndStream(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
        ScalarField3D<T> const &rhoBarField, Dot3D const &offset1,
        TensorField3D<T, Descriptor<T>::d> const &jField, Dot3D const &offset2,
        ScalarField3D<T> const &sigmaField, Dot3D const &offset3, BlockStatistics &stat);
    void boundaryStream(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &bound, Box3D const &domain);

private:
    T rhoBarF, jFsqr;
    Array<T, Descriptor<T>::d> jF;
};

template <typename T, template <typename U> class Descriptor>
class OnLinkExternalRhoJcollideAndStream3D : public BoxProcessingFunctional3D {
public:
    // Block 0: lattice; Block 1: rhoBar; Block 2: j.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks);
    virtual OnLinkExternalRhoJcollideAndStream3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    void collide(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
        ScalarField3D<T> const &rhoBarField, Dot3D const &offset1,
        TensorField3D<T, 3> const &jField, Dot3D const &offset2, BlockStatistics &stat);
    void bulkCollideAndStream(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &domain,
        ScalarField3D<T> const &rhoBarField, Dot3D const &offset1,
        TensorField3D<T, 3> const &jField, Dot3D const &offset2, BlockStatistics &stat);
    void boundaryStream(
        BlockLattice3D<T, Descriptor> &lattice, Box3D const &bound, Box3D const &domain);
};

template <typename T, template <typename U> class Descriptor>
class MaskedCollide3D : public BoxProcessingFunctional3D_LS<T, Descriptor, int> {
public:
    MaskedCollide3D(Dynamics<T, Descriptor> const &dynamics_, int flag_);
    MaskedCollide3D(MaskedCollide3D<T, Descriptor> const &rhs);
    MaskedCollide3D<T, Descriptor> &operator=(MaskedCollide3D<T, Descriptor> const &rhs);
    void swap(MaskedCollide3D<T, Descriptor> &rhs);
    virtual void process(
        Box3D domain, BlockLattice3D<T, Descriptor> &lattice, ScalarField3D<int> &mask);
    virtual MaskedCollide3D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    Dynamics<T, Descriptor> *dynamics;
    int flag;
};

template <typename T, template <typename U> class Descriptor>
void maskedCollide(
    MultiBlockLattice3D<T, Descriptor> &lattice, Dynamics<T, Descriptor> const &dynamics,
    MultiScalarField3D<int> &mask, int flag, Box3D domain);

template <typename T, template <typename U> class Descriptor>
void maskedCollide(
    MultiBlockLattice3D<T, Descriptor> &lattice, Dynamics<T, Descriptor> const &dynamics,
    MultiScalarField3D<int> &mask, int flag);

}  // namespace plb

#endif  // DYNAMICS_PROCESSOR_3D_H
