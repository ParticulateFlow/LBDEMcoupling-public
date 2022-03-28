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
#ifndef FINITE_DIFFERENCE_FUNCTIONAL_2D_HH
#define FINITE_DIFFERENCE_FUNCTIONAL_2D_HH

#include <cmath>

#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "finiteDifference/fdFunctional2D.h"
#include "finiteDifference/fdStencils1D.h"

namespace plb {

/* ******** BoxLaplacianFunctional2D *********************************** */

template <typename T>
void BoxLaplacianFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &laplacian)
{
    Dot2D offset = computeRelativeDisplacement(value, laplacian);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T laplacianValue = -(T)4 * value.get(iX, iY)
                               + (value.get(iX - 1, iY) + value.get(iX + 1, iY)
                                  + value.get(iX, iY - 1) + value.get(iX, iY + 1));
            laplacian.get(iX + offset.x, iY + offset.y) = laplacianValue;
        }
    }
}

template <typename T>
BoxLaplacianFunctional2D<T> *BoxLaplacianFunctional2D<T>::clone() const
{
    return new BoxLaplacianFunctional2D<T>(*this);
}

template <typename T>
void BoxLaplacianFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* ******** BoxXderivativeFunctional2D *********************************** */

template <typename T>
void BoxXderivativeFunctional2D<T>::processBulk(
    Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &derivative)
{
    Dot2D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            derivative.get(iX + offset.x, iY + offset.y) =
                fd::ctl_diff(value.get(iX + 1, iY), value.get(iX - 1, iY));
        }
    }
}

template <typename T>
void BoxXderivativeFunctional2D<T>::processEdge(
    int direction, int orientation, Box2D domain, ScalarField2D<T> &value,
    ScalarField2D<T> &derivative)
{
    Dot2D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            if (direction != 0) {
                derivative.get(iX + offset.x, iY + offset.y) =
                    fd::ctl_diff(value.get(iX + 1, iY), value.get(iX - 1, iY));
            } else {
                if (orientation == 1) {
                    derivative.get(iX + offset.x, iY + offset.y) =
                        fd::o1_fwd_diff(value.get(iX - 1, iY), value.get(iX, iY));
                } else {
                    derivative.get(iX + offset.x, iY + offset.y) =
                        fd::o1_fwd_diff(value.get(iX, iY), value.get(iX + 1, iY));
                }
            }
        }
    }
}

template <typename T>
void BoxXderivativeFunctional2D<T>::processCorner(
    int normalX, int normalY, Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &derivative)
{
    Dot2D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            if (normalX == 1) {
                derivative.get(iX + offset.x, iY + offset.y) =
                    fd::o1_fwd_diff(value.get(iX - 1, iY), value.get(iX, iY));
            } else {
                derivative.get(iX + offset.x, iY + offset.y) =
                    fd::o1_fwd_diff(value.get(iX, iY), value.get(iX + 1, iY));
            }
        }
    }
}

template <typename T>
BoxXderivativeFunctional2D<T> *BoxXderivativeFunctional2D<T>::clone() const
{
    return new BoxXderivativeFunctional2D<T>(*this);
}

template <typename T>
void BoxXderivativeFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxXderivativeFunctional2D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxYderivativeFunctional2D *********************************** */

template <typename T>
void BoxYderivativeFunctional2D<T>::processBulk(
    Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &derivative)
{
    Dot2D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            derivative.get(iX + offset.x, iY + offset.y) =
                fd::ctl_diff(value.get(iX, iY + 1), value.get(iX, iY - 1));
        }
    }
}

template <typename T>
void BoxYderivativeFunctional2D<T>::processEdge(
    int direction, int orientation, Box2D domain, ScalarField2D<T> &value,
    ScalarField2D<T> &derivative)
{
    Dot2D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            if (direction != 1) {
                derivative.get(iX + offset.x, iY + offset.y) =
                    fd::ctl_diff(value.get(iX, iY + 1), value.get(iX, iY - 1));
            } else {
                if (orientation == 1) {
                    derivative.get(iX + offset.x, iY + offset.y) =
                        fd::o1_fwd_diff(value.get(iX, iY - 1), value.get(iX, iY));
                } else {
                    derivative.get(iX + offset.x, iY + offset.y) =
                        fd::o1_fwd_diff(value.get(iX, iY), value.get(iX, iY + 1));
                }
            }
        }
    }
}

template <typename T>
void BoxYderivativeFunctional2D<T>::processCorner(
    int normalX, int normalY, Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &derivative)
{
    Dot2D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            if (normalY == 1) {
                derivative.get(iX + offset.x, iY + offset.y) =
                    fd::o1_fwd_diff(value.get(iX, iY - 1), value.get(iX, iY));
            } else {
                derivative.get(iX + offset.x, iY + offset.y) =
                    fd::o1_fwd_diff(value.get(iX, iY), value.get(iX, iY + 1));
            }
        }
    }
}

template <typename T>
BoxYderivativeFunctional2D<T> *BoxYderivativeFunctional2D<T>::clone() const
{
    return new BoxYderivativeFunctional2D<T>(*this);
}

template <typename T>
void BoxYderivativeFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxYderivativeFunctional2D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxGradientNormFunctional2D *********************************** */

template <typename T>
void BoxGradientNormFunctional2D<T>::processBulk(
    Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &derivative)
{
    Dot2D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T xDeriv = fdDataField::bulkXderiv(value, iX, iY);
            T yDeriv = fdDataField::bulkYderiv(value, iX, iY);
            T gradientNorm = std::sqrt(util::sqr(xDeriv) + util::sqr(yDeriv));
            derivative.get(iX + offset.x, iY + offset.y) = gradientNorm;
        }
    }
}

template <typename T>
void BoxGradientNormFunctional2D<T>::processEdge(
    int direction, int orientation, Box2D domain, ScalarField2D<T> &value,
    ScalarField2D<T> &derivative)
{
    Dot2D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T xDeriv = fdDataField::edgeXderiv(value, direction, orientation, iX, iY);
            T yDeriv = fdDataField::edgeYderiv(value, direction, orientation, iX, iY);
            T gradientNorm = std::sqrt(util::sqr(xDeriv) + util::sqr(yDeriv));
            derivative.get(iX + offset.x, iY + offset.y) = gradientNorm;
        }
    }
}

template <typename T>
void BoxGradientNormFunctional2D<T>::processCorner(
    int normalX, int normalY, Box2D domain, ScalarField2D<T> &value, ScalarField2D<T> &derivative)
{
    Dot2D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T xDeriv = fdDataField::cornerXderiv(value, normalX, normalY, iX, iY);
            T yDeriv = fdDataField::cornerYderiv(value, normalX, normalY, iX, iY);
            T gradientNorm = std::sqrt(util::sqr(xDeriv) + util::sqr(yDeriv));
            derivative.get(iX + offset.x, iY + offset.y) = gradientNorm;
        }
    }
}

template <typename T>
BoxGradientNormFunctional2D<T> *BoxGradientNormFunctional2D<T>::clone() const
{
    return new BoxGradientNormFunctional2D<T>(*this);
}

template <typename T>
void BoxGradientNormFunctional2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxGradientNormFunctional2D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* *************** BoxPoissonResidueFunctional2D ********************* */

template <typename T>
BoxPoissonResidueFunctional2D<T>::BoxPoissonResidueFunctional2D() :
    maxResidueId(this->getStatistics().subscribeMax())
{ }

template <typename T>
void BoxPoissonResidueFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &pressure, ScalarField2D<T> &rhs)
{
    Dot2D offset = computeRelativeDisplacement(pressure, rhs);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sumPressure = pressure.get(iX + 1, iY) + pressure.get(iX, iY + 1)
                            + pressure.get(iX - 1, iY) + pressure.get(iX, iY - 1);
            T residue = std::fabs(
                sumPressure - (T)4 * pressure.get(iX, iY) + rhs.get(iX + offset.x, iY + offset.y));
            statistics.gatherMax(maxResidueId, residue);
        }
    }
}

template <typename T>
BoxPoissonResidueFunctional2D<T> *BoxPoissonResidueFunctional2D<T>::clone() const
{
    return new BoxPoissonResidueFunctional2D<T>(*this);
}

template <typename T>
T BoxPoissonResidueFunctional2D<T>::getMaxResidue() const
{
    return this->getStatistics().getMax(maxResidueId);
}

/* ******** BoxPoissonIteration2D ************************************* */

template <typename T>
BoxPoissonIteration2D<T>::BoxPoissonIteration2D(T beta_) : beta(beta_)
{ }

template <typename T>
void BoxPoissonIteration2D<T>::processBulk(
    Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> &new_u_h = *scalarFields[1];
    ScalarField2D<T> const &rhs = *scalarFields[2];
    Dot2D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sumPressure = u_h.get(iX + 1, iY) + u_h.get(iX, iY + 1) + u_h.get(iX - 1, iY)
                            + u_h.get(iX, iY - 1);

            new_u_h.get(iX + offset1.x, iY + offset1.y) =
                ((T)1 - beta) * u_h.get(iX, iY)
                + (beta / (T)4) * (sumPressure + rhs.get(iX + offset2.x, iY + offset2.y));
        }
    }
}

template <typename T>
void BoxPoissonIteration2D<T>::processEdge(
    int direction, int orientation, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> &new_u_h = *scalarFields[1];
    ScalarField2D<T> const &rhs = *scalarFields[2];
    Dot2D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sumPressure = T();
            if (direction == 0 && orientation == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX + 1, iY);

            if (direction == 1 && orientation == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY + 1);

            if (direction == 0 && orientation == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX - 1, iY);

            if (direction == 1 && orientation == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY - 1);

            new_u_h.get(iX + offset1.x, iY + offset1.y) =
                ((T)1 - beta) * u_h.get(iX, iY)
                + (beta / (T)4) * (sumPressure + rhs.get(iX + offset2.x, iY + offset2.y));
        }
    }
}

template <typename T>
void BoxPoissonIteration2D<T>::processCorner(
    int normalX, int normalY, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> &new_u_h = *scalarFields[1];
    ScalarField2D<T> const &rhs = *scalarFields[2];
    Dot2D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sumPressure = T();
            if (normalX == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX + 1, iY);

            if (normalY == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY + 1);

            if (normalX == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX - 1, iY);

            if (normalY == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY - 1);

            new_u_h.get(iX + offset1.x, iY + offset1.y) =
                ((T)1 - beta) * u_h.get(iX, iY)
                + (beta / (T)4) * (sumPressure + rhs.get(iX + offset2.x, iY + offset2.y));
        }
    }
}

template <typename T>
BoxPoissonIteration2D<T> *BoxPoissonIteration2D<T>::clone() const
{
    return new BoxPoissonIteration2D<T>(*this);
}

template <typename T>
void BoxPoissonIteration2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT BoxPoissonIteration2D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* *************** One Jacobi iteration ************* */
template <typename T>
void JacobiIteration2D<T>::processBulk(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> &new_u_h = *scalarFields[1];
    ScalarField2D<T> const &rhs = *scalarFields[2];
    Dot2D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sum_u_h = u_h.get(iX + 1, iY) + u_h.get(iX, iY + 1) + u_h.get(iX - 1, iY)
                        + u_h.get(iX, iY - 1);

            new_u_h.get(iX + offset1.x, iY + offset1.y) =
                (1.0 / 4.0) * (sum_u_h - rhs.get(iX + offset2.x, iY + offset2.y));
        }
    }
}

template <typename T>
void JacobiIteration2D<T>::processEdge(
    int direction, int orientation, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> &new_u_h = *scalarFields[1];
    ScalarField2D<T> const &rhs = *scalarFields[2];
    Dot2D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sum = T();
            if (direction == 0 && orientation == 1)
                sum += u_h.get(iX, iY);
            else
                sum += u_h.get(iX + 1, iY);

            if (direction == 1 && orientation == 1)
                sum += u_h.get(iX, iY);
            else
                sum += u_h.get(iX, iY + 1);

            if (direction == 0 && orientation == -1)
                sum += u_h.get(iX, iY);
            else
                sum += u_h.get(iX - 1, iY);

            if (direction == 1 && orientation == -1)
                sum += u_h.get(iX, iY);
            else
                sum += u_h.get(iX, iY - 1);

            new_u_h.get(iX + offset1.x, iY + offset1.y) =
                (1.0 / 4.0) * (sum - rhs.get(iX + offset2.x, iY + offset2.y));
        }
    }
}

template <typename T>
void JacobiIteration2D<T>::processCorner(
    int normalX, int normalY, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> &new_u_h = *scalarFields[1];
    ScalarField2D<T> const &rhs = *scalarFields[2];
    Dot2D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sumPressure = T();
            if (normalX == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX + 1, iY);

            if (normalY == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY + 1);

            if (normalX == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX - 1, iY);

            if (normalY == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY - 1);

            new_u_h.get(iX + offset1.x, iY + offset1.y) =
                (1.0 / 4.0) * (sumPressure - rhs.get(iX + offset2.x, iY + offset2.y));
        }
    }
}

template <typename T>
JacobiIteration2D<T> *JacobiIteration2D<T>::clone() const
{
    return new JacobiIteration2D<T>();
}

template <typename T>
void JacobiIteration2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT JacobiIteration2D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* *************** GaussSeidelIteration2D  ****************** */
template <typename T>
void GaussSeidelIteration2D<T>::processBulk(
    Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> const &jacobi_u_h = *scalarFields[1];
    ScalarField2D<T> &new_u_h = *scalarFields[2];
    ScalarField2D<T> const &rhs = *scalarFields[3];
    Dot2D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sumPressure = u_h.get(iX + 1, iY) +         // t
                            u_h.get(iX, iY + 1) +         // t
                            jacobi_u_h.get(iX - 1, iY) +  // t+1
                            jacobi_u_h.get(iX, iY - 1);   // t+1

            new_u_h.get(iX + offset1.x, iY + offset1.y) =
                (1 / (T)4) * (sumPressure - rhs.get(iX + offset2.x, iY + offset2.y));
        }
    }
}

template <typename T>
void GaussSeidelIteration2D<T>::processEdge(
    int direction, int orientation, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> const &jacobi_u_h = *scalarFields[1];
    ScalarField2D<T> &new_u_h = *scalarFields[2];
    ScalarField2D<T> const &rhs = *scalarFields[3];
    Dot2D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sumPressure = T();
            if (direction == 0 && orientation == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX + 1, iY);

            if (direction == 1 && orientation == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY + 1);

            if (direction == 0 && orientation == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += jacobi_u_h.get(iX - 1, iY);

            if (direction == 1 && orientation == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += jacobi_u_h.get(iX, iY - 1);

            new_u_h.get(iX + offset1.x, iY + offset1.y) =
                (1. / 4.) * (sumPressure - rhs.get(iX + offset2.x, iY + offset2.y));
        }
    }
}

template <typename T>
void GaussSeidelIteration2D<T>::processCorner(
    int normalX, int normalY, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> const &jacobi_u_h = *scalarFields[1];
    ScalarField2D<T> &new_u_h = *scalarFields[2];
    ScalarField2D<T> const &rhs = *scalarFields[3];
    Dot2D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sumPressure = T();
            if (normalX == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX + 1, iY);

            if (normalY == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY + 1);

            if (normalX == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += jacobi_u_h.get(iX - 1, iY);

            if (normalY == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += jacobi_u_h.get(iX, iY - 1);

            new_u_h.get(iX + offset1.x, iY + offset1.y) =
                (1. / 4.) * (sumPressure - rhs.get(iX + offset2.x, iY + offset2.y));
        }
    }
}

template <typename T>
GaussSeidelIteration2D<T> *GaussSeidelIteration2D<T>::clone() const
{
    return new GaussSeidelIteration2D<T>();
}

template <typename T>
void GaussSeidelIteration2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT GaussSeidelIteration2D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* *************** Gauss-Seidel defect (d_h)   ****************** */
template <typename T>
void GaussSeidelDefect2D<T>::processBulk(Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> &d_h = *scalarFields[1];
    ScalarField2D<T> const &rhs = *scalarFields[2];

    Dot2D offset1 = computeRelativeDisplacement(u_h, d_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            // Apply laplacien operator to u_h
            T discreteLaplacien_u_h = u_h.get(iX + 1, iY) + u_h.get(iX, iY + 1)
                                      + u_h.get(iX - 1, iY) + u_h.get(iX, iY - 1)
                                      - 4. * u_h.get(iX, iY);
            // compute the defect for the (iX,iY) point
            d_h.get(iX + offset1.x, iY + offset1.y) =
                discreteLaplacien_u_h - rhs.get(iX + offset2.x, iY + offset2.y);
        }
    }
}

template <typename T>
void GaussSeidelDefect2D<T>::processEdge(
    int direction, int orientation, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> &d_h = *scalarFields[1];
    ScalarField2D<T> const &rhs = *scalarFields[2];

    Dot2D offset1 = computeRelativeDisplacement(u_h, d_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sumPressure = T();
            if (direction == 0 && orientation == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX + 1, iY);

            if (direction == 1 && orientation == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY + 1);

            if (direction == 0 && orientation == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX - 1, iY);

            if (direction == 1 && orientation == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY - 1);

            d_h.get(iX + offset1.x, iY + offset1.y) =
                sumPressure - 4.0 * u_h.get(iX, iY) - rhs.get(iX + offset2.x, iY + offset2.y);
        }
    }
}

template <typename T>
void GaussSeidelDefect2D<T>::processCorner(
    int normalX, int normalY, Box2D domain, std::vector<ScalarField2D<T> *> scalarFields)
{
    ScalarField2D<T> const &u_h = *scalarFields[0];
    ScalarField2D<T> &d_h = *scalarFields[1];
    ScalarField2D<T> const &rhs = *scalarFields[2];

    Dot2D offset1 = computeRelativeDisplacement(u_h, d_h);
    Dot2D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T sumPressure = T();
            if (normalX == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX + 1, iY);

            if (normalY == 1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY + 1);

            if (normalX == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX - 1, iY);

            if (normalY == -1)
                sumPressure += u_h.get(iX, iY);
            else
                sumPressure += u_h.get(iX, iY - 1);

            d_h.get(iX + offset1.x, iY + offset1.y) =
                sumPressure - 4.0 * u_h.get(iX, iY) - rhs.get(iX + offset2.x, iY + offset2.y);
        }
    }
}

template <typename T>
GaussSeidelDefect2D<T> *GaussSeidelDefect2D<T>::clone() const
{
    return new GaussSeidelDefect2D();
}

template <typename T>
void GaussSeidelDefect2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT GaussSeidelDefect2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* *************** Max of Gauss-Seidel defect (max(d_h))   ****************** */
template <typename T>
GaussSeidelMaxDefectFunctional2D<T>::GaussSeidelMaxDefectFunctional2D() :
    maxResidueId(this->getStatistics().subscribeMax())
{ }

template <typename T>
void GaussSeidelMaxDefectFunctional2D<T>::process(
    Box2D domain, ScalarField2D<T> &u_h, ScalarField2D<T> &rhs)
{
    Dot2D offset = computeRelativeDisplacement(u_h, rhs);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            T discreteLaplacien_u_h = u_h.get(iX + 1, iY) + u_h.get(iX, iY + 1)
                                      + u_h.get(iX - 1, iY) + u_h.get(iX, iY - 1)
                                      - 4. * u_h.get(iX, iY);
            // compute the absolute value of the residue
            T residue = std::fabs(discreteLaplacien_u_h - rhs.get(iX + offset.x, iY + offset.y));
            // save it in the internal statistics
            statistics.gatherMax(maxResidueId, residue);
        }
    }
}

template <typename T>
GaussSeidelMaxDefectFunctional2D<T> *GaussSeidelMaxDefectFunctional2D<T>::clone() const
{
    return new GaussSeidelMaxDefectFunctional2D<T>(*this);
}

template <typename T>
T GaussSeidelMaxDefectFunctional2D<T>::getMaxResidual() const
{
    return this->getStatistics().getMax(maxResidueId);
}

}  // namespace plb

#endif  // FINITE_DIFFERENCE_FUNCTIONAL_2D_HH
