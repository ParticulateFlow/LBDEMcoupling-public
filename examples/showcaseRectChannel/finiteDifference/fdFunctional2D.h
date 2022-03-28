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
#ifndef FINITE_DIFFERENCE_FUNCTIONAL_2D_H
#define FINITE_DIFFERENCE_FUNCTIONAL_2D_H

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "core/globalDefs.h"

namespace plb {

/* *************** Central finite-difference schemes ***************** */

template <typename T>
class BoxLaplacianFunctional2D : public BoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void process(Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &laplacian);
    virtual BoxLaplacianFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
};

template <typename T>
class BoxXderivativeFunctional2D : public BoundedBoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void processBulk(Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &derivative);
    virtual void processEdge(
        int direction, int orientation, Box2D domain, ScalarField2D<T> &value,
        ScalarField2D<T> &derivative);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, ScalarField2D<T> &value,
        ScalarField2D<T> &derivative);
    virtual BoxXderivativeFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class BoxYderivativeFunctional2D : public BoundedBoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void processBulk(Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &derivative);
    virtual void processEdge(
        int direction, int orientation, Box2D domain, ScalarField2D<T> &value,
        ScalarField2D<T> &derivative);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, ScalarField2D<T> &value,
        ScalarField2D<T> &derivative);
    virtual BoxYderivativeFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class BoxGradientNormFunctional2D : public BoundedBoxProcessingFunctional2D_SS<T, T> {
public:
    virtual void processBulk(Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &grNorm);
    virtual void processEdge(
        int direction, int orientation, Box2D domain, ScalarField2D<T> &value,
        ScalarField2D<T> &grNorm);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &grNorm);
    virtual BoxGradientNormFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** SOR iterations to solve a Poisson equation ******** */

template <typename T>
class BoxPoissonIteration2D : public BoundedScalarFieldBoxProcessingFunctional2D<T> {
public:
    BoxPoissonIteration2D(T beta_);
    virtual void processBulk(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual void processEdge(
        int direction, int orientation, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual BoxPoissonIteration2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;

private:
    T beta;  //< Relaxation parameter
};

/* *************** One Jacobi iteration ************* */
template <typename T>
class JacobiIteration2D : public BoundedScalarFieldBoxProcessingFunctional2D<T> {
public:
    virtual void processBulk(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual void processEdge(
        int direction, int orientation, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual JacobiIteration2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** Gauss-Seidel iterative schema step  ****************** */
template <typename T>
class GaussSeidelIteration2D : public BoundedScalarFieldBoxProcessingFunctional2D<T> {
public:
    virtual void processBulk(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual void processEdge(
        int direction, int orientation, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual GaussSeidelIteration2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** Gauss-Seidel defect (d_h)   ****************** */
/// d_h = discrete_laplacien(u_h) - rhs
/// this defect is important for the multigrid methods
template <typename T>
class GaussSeidelDefect2D : public BoundedScalarFieldBoxProcessingFunctional2D<T> {
public:
    virtual void processBulk(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual void processEdge(
        int direction, int orientation, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual void processCorner(
        int normalX, int normalY, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields);
    virtual GaussSeidelDefect2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

/* *************** Max of Gauss-Seidel defect (max(d_h))   ****************** */
/// if the defect does not need to be computed, we only save the max value in order
/// to analyze the convergence of simple Gauss-Seidel
template <typename T>
class GaussSeidelMaxDefectFunctional2D : public ReductiveBoxProcessingFunctional2D_SS<T, T> {
public:
    GaussSeidelMaxDefectFunctional2D();
    virtual void process(Box2D domain, ScalarField2D<T> &u_h, ScalarField2D<T> &rhs);
    virtual GaussSeidelMaxDefectFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getMaxResidual() const;

private:
    plint maxResidueId;
};

template <typename T>
class BoxPoissonResidueFunctional2D : public ReductiveBoxProcessingFunctional2D_SS<T, T> {
public:
    BoxPoissonResidueFunctional2D();
    virtual void process(Box2D domain, ScalarField2D<T> &pressure, ScalarField2D<T> &rhs);
    virtual BoxPoissonResidueFunctional2D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getMaxResidue() const;

private:
    plint maxResidueId;
};

}  // namespace plb

#endif  // FINITE_DIFFERENCE_FUNCTIONAL_2D_H
