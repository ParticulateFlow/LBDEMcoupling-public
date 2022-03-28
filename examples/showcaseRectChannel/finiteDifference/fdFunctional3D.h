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
#ifndef FINITE_DIFFERENCE_FUNCTIONAL_3D_H
#define FINITE_DIFFERENCE_FUNCTIONAL_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "core/globalDefs.h"

namespace plb {

/* *************** Central finite-difference schemes ***************** */

template <typename T>
class BoxXderivativeFunctional3D : public BoundedBoxProcessingFunctional3D_SS<T, T> {
public:
    virtual void processBulk(Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative);
    virtual void processPlane(
        int direction, int orientation, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &derivative);
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &derivative);
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &derivative);
    virtual BoxXderivativeFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class BoxYderivativeFunctional3D : public BoundedBoxProcessingFunctional3D_SS<T, T> {
public:
    virtual void processBulk(Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative);
    virtual void processPlane(
        int direction, int orientation, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &derivative);
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &derivative);
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &derivative);
    virtual BoxYderivativeFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class BoxZderivativeFunctional3D : public BoundedBoxProcessingFunctional3D_SS<T, T> {
public:
    virtual void processBulk(Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative);
    virtual void processPlane(
        int direction, int orientation, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &derivative);
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &derivative);
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &derivative);
    virtual BoxZderivativeFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class BoxGradientNormFunctional3D : public BoundedBoxProcessingFunctional3D_SS<T, T> {
public:
    virtual void processBulk(Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &grNorm);
    virtual void processPlane(
        int direction, int orientation, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &grNorm);
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &grNorm);
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &value,
        ScalarField3D<T> &grNorm);
    virtual BoxGradientNormFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** SOR iterations to solve a Poisson equation ******** */

template <typename T>
class BoxPoissonIteration3D : public BoundedScalarFieldBoxProcessingFunctional3D<T> {
public:
    BoxPoissonIteration3D(T beta_);
    virtual void processBulk(Box3D domain, std::vector<ScalarField3D<T> *> scalarFields);
    virtual void processPlane(
        int direction, int orientation, Box3D domain, std::vector<ScalarField3D<T> *> scalarFields);
    virtual void processEdge(
        int plane, int normal1, int normal2, Box3D domain,
        std::vector<ScalarField3D<T> *> scalarFields);
    virtual void processCorner(
        int normalX, int normalY, int normalZ, Box3D domain,
        std::vector<ScalarField3D<T> *> scalarFields);
    virtual BoxPoissonIteration3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T beta;  //< Relaxation parameter
};

template <typename T>
class BoxPoissonResidueFunctional3D : public ReductiveBoxProcessingFunctional3D_SS<T, T> {
public:
    BoxPoissonResidueFunctional3D();
    virtual void process(Box3D domain, ScalarField3D<T> &pressure, ScalarField3D<T> &rhs);
    virtual BoxPoissonResidueFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getMaxResidue() const;

private:
    plint maxResidueId;
};

// ========================================================================= //
// PERIODIC VERSIONS OF THE DERIVATIVES AND POISSON SCHEMES //
// ========================================================================= //

/* *************** Central finite-difference schemes ***************** */

template <typename T>
class BoxXperiodicDerivativeFunctional3D : public BoxProcessingFunctional3D_SS<T, T> {
public:
    virtual void process(Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative);
    virtual BoxXperiodicDerivativeFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class BoxYperiodicDerivativeFunctional3D : public BoxProcessingFunctional3D_SS<T, T> {
public:
    virtual void process(Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative);
    virtual BoxYperiodicDerivativeFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class BoxZperiodicDerivativeFunctional3D : public BoxProcessingFunctional3D_SS<T, T> {
public:
    virtual void process(Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative);
    virtual BoxZperiodicDerivativeFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class BoxPeriodicGradientNormFunctional3D : public BoxProcessingFunctional3D_SS<T, T> {
public:
    virtual void process(Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &grNorm);
    virtual BoxPeriodicGradientNormFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class BoxPeriodicGradientFunctional3D : public BoxProcessingFunctional3D_ST<T, T, 3> {
public:
    virtual void process(Box3D domain, ScalarField3D<T> &value, TensorField3D<T, 3> &derivative);
    virtual BoxPeriodicGradientFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** SOR iterations to solve a Poisson equation ******** */

template <typename T>
class BoxPeriodicPoissonIteration3D : public ScalarFieldBoxProcessingFunctional3D<T> {
public:
    BoxPeriodicPoissonIteration3D(T beta_);
    virtual void process(Box3D domain, std::vector<ScalarField3D<T> *> scalarFields);
    virtual BoxPeriodicPoissonIteration3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T beta;  //< Relaxation parameter
};

}  // namespace plb

#endif  // FINITE_DIFFERENCE_FUNCTIONAL_3D_H
