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
#ifndef FINITE_DIFFERENCE_FUNCTIONAL_3D_HH
#define FINITE_DIFFERENCE_FUNCTIONAL_3D_HH

#include <cmath>

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "finiteDifference/fdFunctional3D.h"
#include "finiteDifference/fdStencils1D.h"

namespace plb {

/* ******** BoxXderivativeFunctional3D *********************************** */

template <typename T>
void BoxXderivativeFunctional3D<T>::processBulk(
    Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::bulkXderiv(value, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxXderivativeFunctional3D<T>::processPlane(
    int direction, int orientation, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::planeXderiv(value, direction, orientation, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxXderivativeFunctional3D<T>::processEdge(
    int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::edgeXderiv(value, plane, normal1, normal2, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxXderivativeFunctional3D<T>::processCorner(
    int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::cornerXderiv(value, normalX, normalY, normalZ, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
BoxXderivativeFunctional3D<T> *BoxXderivativeFunctional3D<T>::clone() const
{
    return new BoxXderivativeFunctional3D<T>(*this);
}

template <typename T>
void BoxXderivativeFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxXderivativeFunctional3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxYderivativeFunctional3D *********************************** */

template <typename T>
void BoxYderivativeFunctional3D<T>::processBulk(
    Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::bulkYderiv(value, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxYderivativeFunctional3D<T>::processPlane(
    int direction, int orientation, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::planeYderiv(value, direction, orientation, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxYderivativeFunctional3D<T>::processEdge(
    int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::edgeYderiv(value, plane, normal1, normal2, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxYderivativeFunctional3D<T>::processCorner(
    int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::cornerYderiv(value, normalX, normalY, normalZ, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
BoxYderivativeFunctional3D<T> *BoxYderivativeFunctional3D<T>::clone() const
{
    return new BoxYderivativeFunctional3D<T>(*this);
}

template <typename T>
void BoxYderivativeFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxYderivativeFunctional3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxZderivativeFunctional3D *********************************** */

template <typename T>
void BoxZderivativeFunctional3D<T>::processBulk(
    Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::bulkZderiv(value, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxZderivativeFunctional3D<T>::processPlane(
    int direction, int orientation, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::planeZderiv(value, direction, orientation, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxZderivativeFunctional3D<T>::processEdge(
    int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::edgeZderiv(value, plane, normal1, normal2, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
void BoxZderivativeFunctional3D<T>::processCorner(
    int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::cornerZderiv(value, normalX, normalY, normalZ, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
BoxZderivativeFunctional3D<T> *BoxZderivativeFunctional3D<T>::clone() const
{
    return new BoxZderivativeFunctional3D<T>(*this);
}

template <typename T>
void BoxZderivativeFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxZderivativeFunctional3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxGradientNormFunctional3D *********************************** */

template <typename T>
void BoxGradientNormFunctional3D<T>::processBulk(
    Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T xDeriv = fdDataField::bulkXderiv(value, iX, iY, iZ);
                T yDeriv = fdDataField::bulkYderiv(value, iX, iY, iZ);
                T zDeriv = fdDataField::bulkZderiv(value, iX, iY, iZ);
                T gradientNorm =
                    std::sqrt(util::sqr(xDeriv) + util::sqr(yDeriv) + util::sqr(zDeriv));
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) = gradientNorm;
            }
        }
    }
}

template <typename T>
void BoxGradientNormFunctional3D<T>::processPlane(
    int direction, int orientation, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T xDeriv = fdDataField::planeXderiv(value, direction, orientation, iX, iY, iZ);
                T yDeriv = fdDataField::planeYderiv(value, direction, orientation, iX, iY, iZ);
                T zDeriv = fdDataField::planeZderiv(value, direction, orientation, iX, iY, iZ);
                T gradientNorm =
                    std::sqrt(util::sqr(xDeriv) + util::sqr(yDeriv) + util::sqr(zDeriv));
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) = gradientNorm;
            }
        }
    }
}

template <typename T>
void BoxGradientNormFunctional3D<T>::processEdge(
    int plane, int normal1, int normal2, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T xDeriv = fdDataField::edgeXderiv(value, plane, normal1, normal2, iX, iY, iZ);
                T yDeriv = fdDataField::edgeYderiv(value, plane, normal1, normal2, iX, iY, iZ);
                T zDeriv = fdDataField::edgeZderiv(value, plane, normal1, normal2, iX, iY, iZ);
                T gradientNorm =
                    std::sqrt(util::sqr(xDeriv) + util::sqr(yDeriv) + util::sqr(zDeriv));
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) = gradientNorm;
            }
        }
    }
}

template <typename T>
void BoxGradientNormFunctional3D<T>::processCorner(
    int normalX, int normalY, int normalZ, Box3D domain, ScalarField3D<T> &value,
    ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T xDeriv = fdDataField::cornerXderiv(value, normalX, normalY, normalZ, iX, iY, iZ);
                T yDeriv = fdDataField::cornerYderiv(value, normalX, normalY, normalZ, iX, iY, iZ);
                T zDeriv = fdDataField::cornerZderiv(value, normalX, normalY, normalZ, iX, iY, iZ);
                T gradientNorm =
                    std::sqrt(util::sqr(xDeriv) + util::sqr(yDeriv) + util::sqr(zDeriv));
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) = gradientNorm;
            }
        }
    }
}

template <typename T>
BoxGradientNormFunctional3D<T> *BoxGradientNormFunctional3D<T>::clone() const
{
    return new BoxGradientNormFunctional3D<T>(*this);
}

template <typename T>
void BoxGradientNormFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxGradientNormFunctional3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxPeriodicGradientFunctional3D *********************************** */

template <typename T>
void BoxPeriodicGradientFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &value, TensorField3D<T, 3> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z)[0] =
                    fd::ctl_diff(value.get(iX + 1, iY, iZ), value.get(iX - 1, iY, iZ));
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z)[1] =
                    fd::ctl_diff(value.get(iX, iY + 1, iZ), value.get(iX, iY - 1, iZ));
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z)[2] =
                    fd::ctl_diff(value.get(iX, iY, iZ + 1), value.get(iX, iY, iZ - 1));
            }
        }
    }
}

template <typename T>
BoxPeriodicGradientFunctional3D<T> *BoxPeriodicGradientFunctional3D<T>::clone() const
{
    return new BoxPeriodicGradientFunctional3D<T>(*this);
}

template <typename T>
void BoxPeriodicGradientFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxPeriodicGradientFunctional3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* *************** BoxPoissonResidueFunctional3D ********************* */

template <typename T>
BoxPoissonResidueFunctional3D<T>::BoxPoissonResidueFunctional3D() :
    maxResidueId(this->getStatistics().subscribeMax())
{ }

template <typename T>
void BoxPoissonResidueFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &pressure, ScalarField3D<T> &rhs)
{
    Dot3D offset = computeRelativeDisplacement(pressure, rhs);
    BlockStatistics &statistics = this->getStatistics();
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T sumPressure = pressure.get(iX + 1, iY, iZ) + pressure.get(iX, iY + 1, iZ)
                                + pressure.get(iX, iY, iZ + 1) + pressure.get(iX - 1, iY, iZ)
                                + pressure.get(iX, iY - 1, iZ) + pressure.get(iX, iY, iZ - 1);
                T residue = std::fabs(
                    sumPressure - (T)6 * pressure.get(iX, iY, iZ)
                    + rhs.get(iX + offset.x, iY + offset.y, iZ + offset.z));
                statistics.gatherMax(maxResidueId, residue);
            }
        }
    }
}

template <typename T>
BoxPoissonResidueFunctional3D<T> *BoxPoissonResidueFunctional3D<T>::clone() const
{
    return new BoxPoissonResidueFunctional3D<T>(*this);
}

template <typename T>
T BoxPoissonResidueFunctional3D<T>::getMaxResidue() const
{
    return this->getStatistics().getMax(maxResidueId);
}

/* ******** BoxPoissonIteration3D ************************************* */

template <typename T>
BoxPoissonIteration3D<T>::BoxPoissonIteration3D(T beta_) : beta(beta_)
{ }

template <typename T>
void BoxPoissonIteration3D<T>::processBulk(
    Box3D domain, std::vector<ScalarField3D<T> *> scalarFields)
{
    ScalarField3D<T> const &u_h = *scalarFields[0];
    ScalarField3D<T> &new_u_h = *scalarFields[1];
    ScalarField3D<T> const &rhs = *scalarFields[2];
    Dot3D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot3D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T sumPressure = u_h.get(iX + 1, iY, iZ) + u_h.get(iX, iY + 1, iZ)
                                + u_h.get(iX, iY, iZ + 1) + u_h.get(iX - 1, iY, iZ)
                                + u_h.get(iX, iY - 1, iZ) + u_h.get(iX, iY, iZ - 1);

                new_u_h.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z) =
                    ((T)1 - beta) * u_h.get(iX, iY, iZ)
                    + (beta / (T)6)
                          * (sumPressure + rhs.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z));
            }
        }
    }
}

template <typename T>
void BoxPoissonIteration3D<T>::processPlane(
    int direction, int orientation, Box3D domain, std::vector<ScalarField3D<T> *> scalarFields)
{
    ScalarField3D<T> const &u_h = *scalarFields[0];
    ScalarField3D<T> &new_u_h = *scalarFields[1];
    ScalarField3D<T> const &rhs = *scalarFields[2];
    Dot3D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot3D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T sumPressure = T();
                if (direction == 0 && orientation == 1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX + 1, iY, iZ);

                if (direction == 1 && orientation == 1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX, iY + 1, iZ);

                if (direction == 2 && orientation == 1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX, iY, iZ + 1);

                if (direction == 0 && orientation == -1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX - 1, iY, iZ);

                if (direction == 1 && orientation == -1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX, iY - 1, iZ);

                if (direction == 2 && orientation == -1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX, iY, iZ - 1);

                new_u_h.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z) =
                    ((T)1 - beta) * u_h.get(iX, iY, iZ)
                    + (beta / (T)6)
                          * (sumPressure + rhs.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z));
            }
        }
    }
}

template <typename T>
void BoxPoissonIteration3D<T>::processEdge(
    int plane, int normal1, int normal2, Box3D domain, std::vector<ScalarField3D<T> *> scalarFields)
{
    ScalarField3D<T> const &u_h = *scalarFields[0];
    ScalarField3D<T> &new_u_h = *scalarFields[1];
    ScalarField3D<T> const &rhs = *scalarFields[2];
    Dot3D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot3D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T sumPressure = T();
                if (plane == 0) {
                    sumPressure += u_h.get(iX + 1, iY, iZ);
                    sumPressure += u_h.get(iX - 1, iY, iZ);
                    if (normal1 == 1) {
                        sumPressure += u_h.get(iX, iY - 1, iZ);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }

                    if (normal1 == -1) {
                        sumPressure += u_h.get(iX, iY + 1, iZ);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }

                    if (normal2 == 1) {
                        sumPressure += u_h.get(iX, iY, iZ - 1);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }

                    if (normal2 == -1) {
                        sumPressure += u_h.get(iX, iY, iZ + 1);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }
                } else if (plane == 1) {
                    sumPressure += u_h.get(iX, iY + 1, iZ);
                    sumPressure += u_h.get(iX, iY - 1, iZ);
                    if (normal1 == 1) {
                        sumPressure += u_h.get(iX, iY, iZ - 1);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }

                    if (normal1 == -1) {
                        sumPressure += u_h.get(iX, iY, iZ + 1);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }

                    if (normal2 == 1) {
                        sumPressure += u_h.get(iX - 1, iY, iZ);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }

                    if (normal2 == -1) {
                        sumPressure += u_h.get(iX + 1, iY, iZ);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }
                } else if (plane == 2) {
                    sumPressure += u_h.get(iX, iY, iZ + 1);
                    sumPressure += u_h.get(iX, iY, iZ - 1);
                    if (normal1 == 1) {
                        sumPressure += u_h.get(iX - 1, iY, iZ);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }

                    if (normal1 == -1) {
                        sumPressure += u_h.get(iX + 1, iY, iZ);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }

                    if (normal2 == 1) {
                        sumPressure += u_h.get(iX, iY, iZ - 1);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }

                    if (normal2 == -1) {
                        sumPressure += u_h.get(iX, iY, iZ + 1);
                    } else {
                        sumPressure += u_h.get(iX, iY, iZ);
                    }
                }

                new_u_h.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z) =
                    ((T)1 - beta) * u_h.get(iX, iY, iZ)
                    + (beta / (T)6)
                          * (sumPressure + rhs.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z));
            }
        }
    }
}

template <typename T>
void BoxPoissonIteration3D<T>::processCorner(
    int normalX, int normalY, int normalZ, Box3D domain,
    std::vector<ScalarField3D<T> *> scalarFields)
{
    ScalarField3D<T> const &u_h = *scalarFields[0];
    ScalarField3D<T> &new_u_h = *scalarFields[1];
    ScalarField3D<T> const &rhs = *scalarFields[2];
    Dot3D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot3D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T sumPressure = T();
                if (normalX == 1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX + 1, iY, iZ);

                if (normalY == 1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX, iY + 1, iZ);

                if (normalZ == 1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX, iY, iZ + 1);

                if (normalX == -1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX - 1, iY, iZ);

                if (normalY == -1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX, iY - 1, iZ);

                if (normalZ == -1)
                    sumPressure += u_h.get(iX, iY, iZ);
                else
                    sumPressure += u_h.get(iX, iY, iZ - 1);

                new_u_h.get(iX + offset1.x, iY + offset1.y, iZ + offset1.z) =
                    ((T)1 - beta) * u_h.get(iX, iY, iZ)
                    + (beta / (T)6)
                          * (sumPressure + rhs.get(iX + offset2.x, iY + offset2.y, iZ + offset2.z));
            }
        }
    }
}

template <typename T>
BoxPoissonIteration3D<T> *BoxPoissonIteration3D<T>::clone() const
{
    return new BoxPoissonIteration3D<T>(*this);
}

template <typename T>
void BoxPoissonIteration3D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template <typename T>
BlockDomain::DomainT BoxPoissonIteration3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

// ========================================================================= //
// PERIODIC VERSIONS OF THE DERIVATIVES AND POISSON SCHEMES //
// ========================================================================= //

/* ******** BoxXperiodicDerivativeFunctional3D *********************************** */

template <typename T>
void BoxXperiodicDerivativeFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::bulkXderiv(value, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
BoxXperiodicDerivativeFunctional3D<T> *BoxXperiodicDerivativeFunctional3D<T>::clone() const
{
    return new BoxXperiodicDerivativeFunctional3D<T>(*this);
}

template <typename T>
void BoxXperiodicDerivativeFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxXperiodicDerivativeFunctional3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxYperiodicDerivativeFunctional3D *********************************** */

template <typename T>
void BoxYperiodicDerivativeFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::bulkYderiv(value, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
BoxYperiodicDerivativeFunctional3D<T> *BoxYperiodicDerivativeFunctional3D<T>::clone() const
{
    return new BoxYperiodicDerivativeFunctional3D<T>(*this);
}

template <typename T>
void BoxYperiodicDerivativeFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxYperiodicDerivativeFunctional3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxZperiodicDerivativeFunctional3D *********************************** */

template <typename T>
void BoxZperiodicDerivativeFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) =
                    fdDataField::bulkZderiv(value, iX, iY, iZ);
            }
        }
    }
}

template <typename T>
BoxZperiodicDerivativeFunctional3D<T> *BoxZperiodicDerivativeFunctional3D<T>::clone() const
{
    return new BoxZperiodicDerivativeFunctional3D<T>(*this);
}

template <typename T>
void BoxZperiodicDerivativeFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxZperiodicDerivativeFunctional3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxGradientNormFunctional3D *********************************** */

template <typename T>
void BoxPeriodicGradientNormFunctional3D<T>::process(
    Box3D domain, ScalarField3D<T> &value, ScalarField3D<T> &derivative)
{
    Dot3D offset = computeRelativeDisplacement(value, derivative);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T xDeriv = fdDataField::bulkXderiv(value, iX, iY, iZ);
                T yDeriv = fdDataField::bulkYderiv(value, iX, iY, iZ);
                T zDeriv = fdDataField::bulkZderiv(value, iX, iY, iZ);
                T gradientNorm =
                    std::sqrt(util::sqr(xDeriv) + util::sqr(yDeriv) + util::sqr(zDeriv));
                derivative.get(iX + offset.x, iY + offset.y, iZ + offset.z) = gradientNorm;
            }
        }
    }
}

template <typename T>
BoxPeriodicGradientNormFunctional3D<T> *BoxPeriodicGradientNormFunctional3D<T>::clone() const
{
    return new BoxPeriodicGradientNormFunctional3D<T>(*this);
}

template <typename T>
void BoxPeriodicGradientNormFunctional3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T>
BlockDomain::DomainT BoxPeriodicGradientNormFunctional3D<T>::appliesTo() const
{
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}

/* ******** BoxPeriodicPoissonIteration3D ************************************* */

template <typename T>
BoxPeriodicPoissonIteration3D<T>::BoxPeriodicPoissonIteration3D(T beta_) : beta(beta_)
{ }

template <typename T>
void BoxPeriodicPoissonIteration3D<T>::process(
    Box3D domain, std::vector<ScalarField3D<T> *> scalarFields)
{
    ScalarField3D<T> const &u_h = *scalarFields[0];
    ScalarField3D<T> &new_u_h = *scalarFields[1];
    ScalarField3D<T> const &rhs = *scalarFields[2];
    Dot3D offset1 = computeRelativeDisplacement(u_h, new_u_h);
    Dot3D offset2 = computeRelativeDisplacement(u_h, rhs);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint oX1 = iX + offset1.x;
        plint oX2 = iX + offset2.x;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint oY1 = iY + offset1.y;
            plint oY2 = iY + offset2.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T sumPressure = u_h.get(iX + 1, iY, iZ) + u_h.get(iX, iY + 1, iZ)
                                + u_h.get(iX, iY, iZ + 1) + u_h.get(iX - 1, iY, iZ)
                                + u_h.get(iX, iY - 1, iZ) + u_h.get(iX, iY, iZ - 1);

                new_u_h.get(oX1, oY1, iZ + offset1.z) =
                    ((T)1 - beta) * u_h.get(iX, iY, iZ)
                    + (beta / (T)6) * (sumPressure + rhs.get(oX2, oY2, iZ + offset2.z));
            }
        }
    }
}

template <typename T>
BoxPeriodicPoissonIteration3D<T> *BoxPeriodicPoissonIteration3D<T>::clone() const
{
    return new BoxPeriodicPoissonIteration3D<T>(*this);
}

template <typename T>
void BoxPeriodicPoissonIteration3D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

}  // namespace plb

#endif  // FINITE_DIFFERENCE_FUNCTIONAL_3D_HH
