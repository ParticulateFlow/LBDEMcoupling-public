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

#ifndef FREE_SURFACE_MODEL_3D_HH
#define FREE_SURFACE_MODEL_3D_HH

#include <limits>

#include "atomicBlock/atomicContainerBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/block3D.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiPhysics/freeSurfaceModel3D.h"
#include "multiPhysics/freeSurfaceTemplates.h"

namespace plb {

/* *************** Class FreeSurfaceComputeNormals3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceComputeNormals3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    // Smooth the volume fraction twice. (At the end include also a 1-cell layer around "domain".)
    plint nx = domain.getNx() + 4;
    plint ny = domain.getNy() + 4;
    plint nz = domain.getNz() + 4;
    ScalarField3D<T> smoothVolumeFractionTmp(nx, ny, nz);
    Dot3D ofsSVFT(-domain.x0 + 2, -domain.y0 + 2, -domain.z0 + 2);
    for (plint iX = domain.x0 - 2; iX <= domain.x1 + 2; ++iX) {
        plint i = iX + ofsSVFT.x;
        for (plint iY = domain.y0 - 2; iY <= domain.y1 + 2; ++iY) {
            plint j = iY + ofsSVFT.y;
            for (plint iZ = domain.z0 - 2; iZ <= domain.z1 + 2; ++iZ) {
                plint k = iZ + ofsSVFT.z;
                // smoothVolumeFractionTmp.get(i, j, k) = param.smooth(*param.volumeFractionP(),
                //         param.volumeFractionOffset(), iX, iY, iZ);
                // smoothVolumeFractionTmp.get(i, j, k) =
                //     param.template lbmSmooth<T, Descriptor>(
                //             *param.volumeFractionP(), param.volumeFractionOffset(), iX, iY, iZ);
                smoothVolumeFractionTmp.get(i, j, k) =
                    param.template lbmSmooth<T, descriptors::AdvectionDiffusionD3Q7Descriptor>(
                        *param.volumeFractionP(), param.volumeFractionOffset(), iX, iY, iZ);
            }
        }
    }

    nx = domain.getNx() + 2;
    ny = domain.getNy() + 2;
    nz = domain.getNz() + 2;
    ScalarField3D<T> smoothVolumeFraction(nx, ny, nz);
    Dot3D ofsSVF(-domain.x0 + 1, -domain.y0 + 1, -domain.z0 + 1);
    for (plint iX = domain.x0 - 1; iX <= domain.x1 + 1; ++iX) {
        plint i = iX + ofsSVF.x;
        for (plint iY = domain.y0 - 1; iY <= domain.y1 + 1; ++iY) {
            plint j = iY + ofsSVF.y;
            for (plint iZ = domain.z0 - 1; iZ <= domain.z1 + 1; ++iZ) {
                plint k = iZ + ofsSVF.z;
                // smoothVolumeFraction.get(i, j, k) = param.smooth(smoothVolumeFractionTmp,
                // ofsSVFT, iX, iY, iZ); smoothVolumeFraction.get(i, j, k) =
                //     param.template lbmSmooth<T, Descriptor>(
                //             smoothVolumeFractionTmp, ofsSVFT, iX, iY, iZ);
                smoothVolumeFraction.get(i, j, k) =
                    param.template lbmSmooth<T, descriptors::AdvectionDiffusionD3Q7Descriptor>(
                        smoothVolumeFractionTmp, ofsSVFT, iX, iY, iZ);
            }
        }
    }

    // The outward pointing unit normal is: n = - grad(VOF) / ||grad(VOF)||.

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Array<T, 3> normal((T)0, (T)0, (T)0);

                if (isAnyWall(param.flag(iX, iY, iZ))) {
                    param.setNormal(iX, iY, iZ, normal);
                    continue;
                }

                /*
                int useLB = 1;
                for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if (isAnyWall(param.flag(nextX, nextY, nextZ))) {
                        useLB = 0;
                        break;
                    }
                }

                if (useLB) {
                */
                // Compute the gradient of the smoothed volume fraction "the lattice Boltzmann way".
                normal = param.template lbmComputeGradient<T, Descriptor>(
                    smoothVolumeFraction, ofsSVF, iX, iY, iZ);
                /*
                } else {
                    // Compute the gradient of the smoothed volume fraction with finite differences
                    // excluding the wall cells.
                    plint h = 1;    // This must be 1 because above we included only a 1-cell layer
                around "domain". normal = param.computeGradient(smoothVolumeFraction, ofsSVF, h, iX,
                iY, iZ);
                }
                */

                T nn = norm(normal);
                if (util::isZero(nn)) {
                    normal.resetToZero();
                } else {
                    normal /= -nn;
                }
                param.setNormal(iX, iY, iZ, normal);
            }
        }
    }
}

/* *************** Class FreeSurfaceGeometry3D ******************************** */

template <typename T, template <typename U> class Descriptor>
ScalarField3D<int> *FreeSurfaceGeometry3D<T, Descriptor>::getInterfaceFlags(
    Box3D domain, FreeSurfaceProcessorParam3D<T, Descriptor> &param)
{
    using namespace freeSurfaceFlag;

    // Define a temporary scalar field for local use in this function. This scalar field will
    // contain 1 extra layer of cells around "domain".
    plint nx = domain.x1 - domain.x0 + 1;
    plint ny = domain.y1 - domain.y0 + 1;
    plint nz = domain.z1 - domain.z0 + 1;
    ScalarField3D<int> *tmp = new ScalarField3D<int>(nx + 2, ny + 2, nz + 2, (int)unTagged);
    PLB_ASSERT(tmp);

    // First tag all regular and contact line interface cells. (Loop along 1 envelope cell as well).
    // All interface tags are stored in the temporary storage.
    for (plint iX = domain.x0 - 1; iX <= domain.x1 + 1; ++iX) {
        plint indX = iX - domain.x0 + 1;
        for (plint iY = domain.y0 - 1; iY <= domain.y1 + 1; ++iY) {
            plint indY = iY - domain.y0 + 1;
            for (plint iZ = domain.z0 - 1; iZ <= domain.z1 + 1; ++iZ) {
                plint indZ = iZ - domain.z0 + 1;
                if (param.flag(iX, iY, iZ) != interface) {
                    tmp->get(indX, indY, indZ) = notInterface;
                    continue;
                }

                // Find all wall neighbors and store their indices.
                int numWallNeighbors = 0;
                std::vector<Array<plint, 3> > wallNeighborIndex;
                for (int dx = -1; dx < 2; dx++) {
                    plint i = iX + dx;
                    for (int dy = -1; dy < 2; dy++) {
                        plint j = iY + dy;
                        for (int dz = -1; dz < 2; dz++) {
                            plint k = iZ + dz;
                            if (!(dx == 0 && dy == 0 && dz == 0)) {
                                if (isAnyWall(param.flag(i, j, k))) {
                                    numWallNeighbors++;
                                    wallNeighborIndex.push_back(Array<plint, 3>(i, j, k));
                                }
                            }
                        }
                    }
                }

                if (numWallNeighbors == 0) {
                    tmp->get(indX, indY, indZ) = regular;
                    continue;
                }

                for (int dx = -1; dx < 2; dx++) {
                    plint i = iX + dx;
                    for (int dy = -1; dy < 2; dy++) {
                        plint j = iY + dy;
                        for (int dz = -1; dz < 2; dz++) {
                            plint k = iZ + dz;
                            if (!(dx == 0 && dy == 0 && dz == 0)) {
                                if ((contactAngle > 90.0 && isFullWet(param.flag(i, j, k)))
                                    || (contactAngle <= 90.0 && isEmpty(param.flag(i, j, k))))
                                {
                                    for (int dxx = -1; dxx < 2; dxx++) {
                                        plint ii = i + dxx;
                                        for (int dyy = -1; dyy < 2; dyy++) {
                                            plint jj = j + dyy;
                                            for (int dzz = -1; dzz < 2; dzz++) {
                                                plint kk = k + dzz;
                                                if (!(dxx == 0 && dyy == 0 && dzz == 0)) {
                                                    if (isAnyWall(param.flag(ii, jj, kk))) {
                                                        for (int iWall = 0;
                                                             iWall < numWallNeighbors; iWall++) {
                                                            if (ii == wallNeighborIndex[iWall][0]
                                                                && jj == wallNeighborIndex[iWall][1]
                                                                && kk
                                                                       == wallNeighborIndex[iWall]
                                                                                           [2])
                                                            {
                                                                tmp->get(indX, indY, indZ) =
                                                                    contactLine;
                                                                goto label0;
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
label0:
                continue;
            }
        }
    }

    // Define a scalar field with the interface flags that will be returned from this function.
    ScalarField3D<int> *interfaceFlag = new ScalarField3D<int>(nx, ny, nz, (int)unTagged);
    PLB_ASSERT(interfaceFlag);

    // Now tag all adjacent interface cells and copy all information to the scalar field to be
    // returned. At this point all cells that have the flag "unTagged" are non-contact-line
    // interface cells with wall neighbors, so they are either "regular" or "adjacent" interface
    // cells.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        plint indXtmp = iX - domain.x0 + 1;
        plint indX = iX - domain.x0;
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            plint indYtmp = iY - domain.y0 + 1;
            plint indY = iY - domain.y0;
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                plint indZtmp = iZ - domain.z0 + 1;
                plint indZ = iZ - domain.z0;
                if (tmp->get(indXtmp, indYtmp, indZtmp) != unTagged) {
                    interfaceFlag->get(indX, indY, indZ) = tmp->get(indXtmp, indYtmp, indZtmp);
                } else {
                    int isAdjacent = 0;
                    for (int dx = -1; dx < 2; dx++) {
                        plint i = indXtmp + dx;
                        for (int dy = -1; dy < 2; dy++) {
                            plint j = indYtmp + dy;
                            for (int dz = -1; dz < 2; dz++) {
                                plint k = indZtmp + dz;
                                if (!(dx == 0 && dy == 0 && dz == 0)) {
                                    if (tmp->get(i, j, k) == contactLine) {
                                        isAdjacent = 1;
                                        interfaceFlag->get(indX, indY, indZ) = adjacent;
                                        goto label1;
                                    }
                                }
                            }
                        }
                    }
label1:
                    if (!isAdjacent) {
                        interfaceFlag->get(indX, indY, indZ) = regular;
                    }
                }
            }
        }
    }

    // Check for untagged cells
#ifdef PLB_DEBUG
    for (plint i = 0; i < nx; i++) {
        for (plint j = 0; j < ny; j++) {
            for (plint k = 0; k < nz; k++) {
                PLB_ASSERT(interfaceFlag->get(i, j, k) != unTagged);
            }
        }
    }
#endif

    delete tmp;

    return interfaceFlag;
}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceGeometry3D<T, Descriptor>::computeHeights3D(
    FreeSurfaceProcessorParam3D<T, Descriptor> &param, int integrationDirection, plint iX, plint iY,
    plint iZ, T h[3][3])
{
    using namespace freeSurfaceFlag;

    // Compute the vector parallel to the integration direction.
    Array<int, 3> integrationVector;
    integrationVector[0] = integrationDirection == 0 ? 1 : 0;
    integrationVector[1] = integrationDirection == 1 ? 1 : 0;
    integrationVector[2] = integrationDirection == 2 ? 1 : 0;

    // Compute the vectors tangent to the plane which is normal to the integration vector.
    int iTangentDirection0 = integrationDirection == 0 ? 1 : (integrationDirection == 1) ? 2 : 0;
    int iTangentDirection1 = integrationDirection == 0 ? 2 : (integrationDirection == 1) ? 0 : 1;
    Array<int, 3> tangent0;
    tangent0[0] = iTangentDirection0 == 0 ? 1 : 0;
    tangent0[1] = iTangentDirection0 == 1 ? 1 : 0;
    tangent0[2] = iTangentDirection0 == 2 ? 1 : 0;
    Array<int, 3> tangent1;
    tangent1[0] = iTangentDirection1 == 0 ? 1 : 0;
    tangent1[1] = iTangentDirection1 == 1 ? 1 : 0;
    tangent1[2] = iTangentDirection1 == 2 ? 1 : 0;

    // Calculate the integration stencil width.
    int maxLim = 3;
    for (int d0 = -1; d0 <= 1; d0++) {
        for (int d1 = -1; d1 <= 1; d1++) {
            plint posX = iX + d0 * tangent0[0] + d1 * tangent1[0];
            plint posY = iY + d0 * tangent0[1] + d1 * tangent1[1];
            plint posZ = iZ + d0 * tangent0[2] + d1 * tangent1[2];
            if (isAnyWall(param.flag(posX, posY, posZ))) {
                continue;
            }
            for (int d = 1; d <= maxLim; d++) {
                plint nextX = posX + d * integrationVector[0];
                plint nextY = posY + d * integrationVector[1];
                plint nextZ = posZ + d * integrationVector[2];
                if (isAnyWall(param.flag(nextX, nextY, nextZ))) {
                    maxLim = std::min(maxLim, d - 1);
                    break;
                }
            }
        }
    }

    int minLim = 3;
    for (int d0 = -1; d0 <= 1; d0++) {
        for (int d1 = -1; d1 <= 1; d1++) {
            plint posX = iX + d0 * tangent0[0] + d1 * tangent1[0];
            plint posY = iY + d0 * tangent0[1] + d1 * tangent1[1];
            plint posZ = iZ + d0 * tangent0[2] + d1 * tangent1[2];
            if (isAnyWall(param.flag(posX, posY, posZ))) {
                continue;
            }
            for (int d = 1; d <= minLim; d++) {
                plint nextX = posX - d * integrationVector[0];
                plint nextY = posY - d * integrationVector[1];
                plint nextZ = posZ - d * integrationVector[2];
                if (isAnyWall(param.flag(nextX, nextY, nextZ))) {
                    minLim = std::min(minLim, d - 1);
                    break;
                }
            }
        }
    }

    // Properly initialize heights to -1.
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            h[i][j] = -1.0;
        }
    }

    // Integrate.
    for (int d0 = -1; d0 <= 1; d0++) {
        int i = d0 + 1;
        for (int d1 = -1; d1 <= 1; d1++) {
            int j = d1 + 1;

            plint posX = iX + d0 * tangent0[0] + d1 * tangent1[0];
            plint posY = iY + d0 * tangent0[1] + d1 * tangent1[1];
            plint posZ = iZ + d0 * tangent0[2] + d1 * tangent1[2];
            if (isAnyWall(param.flag(posX, posY, posZ))) {
                continue;
            }
            h[i][j] = 0.0;
            for (int d = -minLim; d <= maxLim; d++) {
                plint nextX = posX + d * integrationVector[0];
                plint nextY = posY + d * integrationVector[1];
                plint nextZ = posZ + d * integrationVector[2];
                h[i][j] += param.volumeFraction(nextX, nextY, nextZ);
            }
        }
    }

    // Extrapolate on walls. (No contact angle algorithm).
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (util::isZero(h[i][j] + (T)1)) {
                h[i][j] = h[1][1];
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceGeometry3D<T, Descriptor>::computeHeights2D(
    FreeSurfaceProcessorParam3D<T, Descriptor> &param, Array<int, 3> &wallTangent0,
    Array<int, 3> &wallTangent1, int integrationDirection, plint iX, plint iY, plint iZ, T h[3])
{
    using namespace freeSurfaceFlag;

    // Compute the vector parallel to the integration direction.
    Array<int, 3> integrationVector;
    integrationVector[0] = integrationDirection == 0 ? wallTangent0[0] : wallTangent1[0];
    integrationVector[1] = integrationDirection == 0 ? wallTangent0[1] : wallTangent1[1];
    integrationVector[2] = integrationDirection == 0 ? wallTangent0[2] : wallTangent1[2];

    // Compute the vector tangent to the line which is normal to the integration vector.
    Array<int, 3> tangent;
    tangent[0] = integrationDirection == 0 ? wallTangent1[0] : wallTangent0[0];
    tangent[1] = integrationDirection == 0 ? wallTangent1[1] : wallTangent0[1];
    tangent[2] = integrationDirection == 0 ? wallTangent1[2] : wallTangent0[2];

    // Calculate the integration stencil width.
    int maxLim = 3;
    for (int d0 = -1; d0 <= 1; d0++) {
        plint posX = iX + d0 * tangent[0];
        plint posY = iY + d0 * tangent[1];
        plint posZ = iZ + d0 * tangent[2];
        if (isAnyWall(param.flag(posX, posY, posZ))) {
            continue;
        }
        for (int d = 1; d <= maxLim; d++) {
            plint nextX = posX + d * integrationVector[0];
            plint nextY = posY + d * integrationVector[1];
            plint nextZ = posZ + d * integrationVector[2];
            if (isAnyWall(param.flag(nextX, nextY, nextZ))) {
                maxLim = std::min(maxLim, d - 1);
                break;
            }
        }
    }

    int minLim = 3;
    for (int d0 = -1; d0 <= 1; d0++) {
        plint posX = iX + d0 * tangent[0];
        plint posY = iY + d0 * tangent[1];
        plint posZ = iZ + d0 * tangent[2];
        if (isAnyWall(param.flag(posX, posY, posZ))) {
            continue;
        }
        for (int d = 1; d <= minLim; d++) {
            plint nextX = posX - d * integrationVector[0];
            plint nextY = posY - d * integrationVector[1];
            plint nextZ = posZ - d * integrationVector[2];
            if (isAnyWall(param.flag(nextX, nextY, nextZ))) {
                minLim = std::min(minLim, d - 1);
                break;
            }
        }
    }

    // Properly initialize heights to -1.
    h[0] = h[1] = h[2] = -1.0;

    // Integrate.
    for (int d0 = -1; d0 <= 1; d0++) {
        int i = d0 + 1;
        plint posX = iX + d0 * tangent[0];
        plint posY = iY + d0 * tangent[1];
        plint posZ = iZ + d0 * tangent[2];
        if (isAnyWall(param.flag(posX, posY, posZ))) {
            continue;
        }
        h[i] = 0.0;
        for (int d = -minLim; d <= maxLim; d++) {
            plint nextX = posX + d * integrationVector[0];
            plint nextY = posY + d * integrationVector[1];
            plint nextZ = posZ + d * integrationVector[2];
            h[i] += param.volumeFraction(nextX, nextY, nextZ);
        }
    }

    // Extrapolate on walls. (No contact angle algorithm).
    for (int i = 0; i < 3; i++) {
        if (util::isZero(h[i] + (T)1)) {
            h[i] = h[1];
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceGeometry3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    static T degToRad = (T)3.14159265358979323844L / (T)180;

    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    Array<T, 3> zeroVector((T)0, (T)0, (T)0);

    if (!useContactAngle) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    if (param.flag(iX, iY, iZ) == interface) {
                        if (param.isBoundary(iX, iY, iZ)) {
                            param.curvature(iX, iY, iZ) = 0.0;
                            param.setNormal(iX, iY, iZ, zeroVector);
                            continue;
                        }
                        // Locally smooth the volume fraction to compute an estimate of the normal.
                        T svfcp = param.smooth(
                            *param.volumeFractionP(), param.volumeFractionOffset(), iX, iY, iZ);
                        T svfx0 = !isAnyWall(param.flag(iX - 1, iY, iZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(),
                                      iX - 1, iY, iZ)
                                                                         : svfcp;
                        T svfx1 = !isAnyWall(param.flag(iX + 1, iY, iZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(),
                                      iX + 1, iY, iZ)
                                                                         : svfcp;
                        T svfy0 = !isAnyWall(param.flag(iX, iY - 1, iZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), iX,
                                      iY - 1, iZ)
                                                                         : svfcp;
                        T svfy1 = !isAnyWall(param.flag(iX, iY + 1, iZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), iX,
                                      iY + 1, iZ)
                                                                         : svfcp;
                        T svfz0 = !isAnyWall(param.flag(iX, iY, iZ - 1)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), iX,
                                      iY, iZ - 1)
                                                                         : svfcp;
                        T svfz1 = !isAnyWall(param.flag(iX, iY, iZ + 1)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), iX,
                                      iY, iZ + 1)
                                                                         : svfcp;

                        // Compute a normalized grad(VF) (inward-pointing normal).
                        Array<T, 3> gradVF;
                        gradVF[0] = 0.5 * (svfx1 - svfx0);
                        gradVF[1] = 0.5 * (svfy1 - svfy0);
                        gradVF[2] = 0.5 * (svfz1 - svfz0);
                        T norm_gradVF = norm(gradVF);
                        if (util::isZero(norm_gradVF)) {
                            param.curvature(iX, iY, iZ) = 0.0;
                            param.setNormal(iX, iY, iZ, zeroVector);
                            continue;
                        }
                        gradVF /= norm_gradVF;

                        T abs0 = std::fabs(gradVF[0]);
                        T abs1 = std::fabs(gradVF[1]);
                        T abs2 = std::fabs(gradVF[2]);

                        int integrationDirection = 2;

                        if (abs0 > abs1) {
                            if (abs0 > abs2) {
                                integrationDirection = 0;
                            }
                            // abs0>abs1 && abs0 <= abs2
                            else
                            {
                                integrationDirection = 2;
                            }
                        }
                        // abs0 <= abs1
                        else
                        {
                            if (abs1 > abs2) {
                                integrationDirection = 1;
                            }
                            // abs0 <= abs1 && abs1 <= abs2
                            else
                            {
                                integrationDirection = 2;
                            }
                        }

                        T h[3][3];
                        computeHeights3D(param, integrationDirection, iX, iY, iZ, h);

                        T dh0 = 0.5 * (h[2][1] - h[0][1]);
                        T dh1 = 0.5 * (h[1][2] - h[1][0]);

                        T dh00 = h[2][1] - 2.0 * h[1][1] + h[0][1];
                        T dh11 = h[1][2] - 2.0 * h[1][1] + h[1][0];

                        T dh01 = 0.25 * (h[2][2] - h[2][0] - h[0][2] + h[0][0]);

                        T value = -(dh00 + dh11 + dh00 * dh1 * dh1 + dh11 * dh0 * dh0
                                    - 2.0 * dh01 * dh0 * dh1)
                                  / std::pow((T)1.0 + dh0 * dh0 + dh1 * dh1, (T)1.5);

                        param.curvature(iX, iY, iZ) = value;

                        T sgn = -gradVF[integrationDirection] < 0.0 ? -1.0 : 1.0;
                        Array<T, 3> normal;
                        if (integrationDirection == 0) {
                            normal = Array<T, 3>(sgn, -dh0, -dh1);
                        } else if (integrationDirection == 1) {
                            normal = Array<T, 3>(-dh1, sgn, -dh0);
                        } else {
                            normal = Array<T, 3>(-dh0, -dh1, sgn);
                        }
                        T norm_normal = norm(normal);
                        if (util::isZero(norm_normal)) {
                            param.setNormal(iX, iY, iZ, zeroVector);
                        } else {
                            param.setNormal(iX, iY, iZ, normal / norm_normal);
                        }
                    } else {
                        param.curvature(iX, iY, iZ) = 0.0;
                        param.setNormal(iX, iY, iZ, zeroVector);
                    }
                }
            }
        }
    } else {  // Use contact angles.
        // First compute the flags.
        ScalarField3D<int> *interfaceFlag = getInterfaceFlags(domain, param);

        /* New contact angle algorithm. This algorithm does not treat properly either the adjacent
         * cells or the free-slip wall cells. */

        // First loop over all the regular and adjacent interface cells and calculate the curvature
        // and the normal vectors. When the appropriate algorithm for the adjacent cells is
        // implemented, they must be removed from these loops.
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX - domain.x0;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY - domain.y0;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint k = iZ - domain.z0;
                    if (interfaceFlag->get(i, j, k) == regular
                        || interfaceFlag->get(i, j, k) == adjacent) {
                        if (param.isBoundary(iX, iY, iZ)) {
                            param.curvature(iX, iY, iZ) = 0.0;
                            param.setNormal(iX, iY, iZ, zeroVector);
                            continue;
                        }
                        // Locally smooth the volume fraction to compute an estimate of the normal.
                        T svfcp = param.smooth(
                            *param.volumeFractionP(), param.volumeFractionOffset(), iX, iY, iZ);
                        T svfx0 = !isAnyWall(param.flag(iX - 1, iY, iZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(),
                                      iX - 1, iY, iZ)
                                                                         : svfcp;
                        T svfx1 = !isAnyWall(param.flag(iX + 1, iY, iZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(),
                                      iX + 1, iY, iZ)
                                                                         : svfcp;
                        T svfy0 = !isAnyWall(param.flag(iX, iY - 1, iZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), iX,
                                      iY - 1, iZ)
                                                                         : svfcp;
                        T svfy1 = !isAnyWall(param.flag(iX, iY + 1, iZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), iX,
                                      iY + 1, iZ)
                                                                         : svfcp;
                        T svfz0 = !isAnyWall(param.flag(iX, iY, iZ - 1)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), iX,
                                      iY, iZ - 1)
                                                                         : svfcp;
                        T svfz1 = !isAnyWall(param.flag(iX, iY, iZ + 1)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), iX,
                                      iY, iZ + 1)
                                                                         : svfcp;

                        // Compute a normalized grad(VF) (inward-pointing normal).
                        Array<T, 3> gradVF;
                        gradVF[0] = 0.5 * (svfx1 - svfx0);
                        gradVF[1] = 0.5 * (svfy1 - svfy0);
                        gradVF[2] = 0.5 * (svfz1 - svfz0);
                        T norm_gradVF = norm(gradVF);
                        if (util::isZero(norm_gradVF)) {
                            param.curvature(iX, iY, iZ) = 0.0;
                            param.setNormal(iX, iY, iZ, zeroVector);
                            continue;
                        }
                        gradVF /= norm_gradVF;

                        T abs0 = std::fabs(gradVF[0]);
                        T abs1 = std::fabs(gradVF[1]);
                        T abs2 = std::fabs(gradVF[2]);

                        int integrationDirection = 2;

                        if (abs0 > abs1) {
                            if (abs0 > abs2) {
                                integrationDirection = 0;
                            }
                            // abs0>abs1 && abs0 <= abs2
                            else
                            {
                                integrationDirection = 2;
                            }
                        }
                        // abs0 <= abs1
                        else
                        {
                            if (abs1 > abs2) {
                                integrationDirection = 1;
                            }
                            // abs0 <= abs1 && abs1 <= abs2
                            else
                            {
                                integrationDirection = 2;
                            }
                        }

                        T h[3][3];
                        computeHeights3D(param, integrationDirection, iX, iY, iZ, h);

                        T dh0 = 0.5 * (h[2][1] - h[0][1]);
                        T dh1 = 0.5 * (h[1][2] - h[1][0]);

                        T dh00 = h[2][1] - 2.0 * h[1][1] + h[0][1];
                        T dh11 = h[1][2] - 2.0 * h[1][1] + h[1][0];

                        T dh01 = 0.25 * (h[2][2] - h[2][0] - h[0][2] + h[0][0]);

                        T value = -(dh00 + dh11 + dh00 * dh1 * dh1 + dh11 * dh0 * dh0
                                    - 2.0 * dh01 * dh0 * dh1)
                                  / std::pow((T)1.0 + dh0 * dh0 + dh1 * dh1, (T)1.5);

                        param.curvature(iX, iY, iZ) = value;

                        T sgn = -gradVF[integrationDirection] < 0.0 ? -1.0 : 1.0;
                        Array<T, 3> normal;
                        if (integrationDirection == 0) {
                            normal = Array<T, 3>(sgn, -dh0, -dh1);
                        } else if (integrationDirection == 1) {
                            normal = Array<T, 3>(-dh1, sgn, -dh0);
                        } else {
                            normal = Array<T, 3>(-dh0, -dh1, sgn);
                        }
                        T norm_normal = norm(normal);
                        if (util::isZero(norm_normal)) {
                            param.setNormal(iX, iY, iZ, zeroVector);
                        } else {
                            param.setNormal(iX, iY, iZ, normal / norm_normal);
                        }
                    } else {
                        param.curvature(iX, iY, iZ) = 0.0;
                        param.setNormal(iX, iY, iZ, zeroVector);
                    }
                }
            }
        }

        // Then loop over all the contact-line interface cells and calculate the curvature and the
        // normal vectors according to the specified contact angle.
        T contactAngleRad = contactAngle * degToRad;
        T tanContactAngle = std::tan(contactAngleRad);

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint i = iX - domain.x0;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint j = iY - domain.y0;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint k = iZ - domain.z0;
                    if (interfaceFlag->get(i, j, k) == contactLine) {
                        if (param.isBoundary(iX, iY, iZ)) {
                            param.curvature(iX, iY, iZ) = 0.0;
                            param.setNormal(iX, iY, iZ, zeroVector);
                            continue;
                        }

                        // First decide where is the wall.
                        int numWallCells = 0;
                        // Computation of the inward-pointing wall normal (not unitary).
                        Array<int, 3> inwardWallNormal(0, 0, 0);
                        for (int dx = -1; dx < 2; dx++) {
                            for (int dy = -1; dy < 2; dy++) {
                                for (int dz = -1; dz < 2; dz++) {
                                    if (isAnyWall(param.flag(iX + dx, iY + dy, iZ + dz))) {
                                        inwardWallNormal += Array<int, 3>(-dx, -dy, -dz);
                                        numWallCells++;
                                    }
                                }
                            }
                        }
                        PLB_ASSERT(numWallCells != 0);
#ifdef PLB_DEBUG
                        int norm2_inwardWallNormal = inwardWallNormal[0] * inwardWallNormal[0]
                                                     + inwardWallNormal[1] * inwardWallNormal[1]
                                                     + inwardWallNormal[2] * inwardWallNormal[2];
#endif
                        PLB_ASSERT(norm2_inwardWallNormal != 0);

                        int iWallNormalDirection;
                        // The inwardWallNormal is aligned with one axis.
                        if (inwardWallNormal[0] != 0 && inwardWallNormal[1] == 0
                            && inwardWallNormal[2] == 0) {
                            iWallNormalDirection = 0;
                        } else if (
                            inwardWallNormal[1] != 0 && inwardWallNormal[2] == 0
                            && inwardWallNormal[0] == 0) {
                            iWallNormalDirection = 1;
                        } else if (
                            inwardWallNormal[2] != 0 && inwardWallNormal[0] == 0
                            && inwardWallNormal[1] == 0) {
                            iWallNormalDirection = 2;
                        } else {
                            // The inwardWallNormal is not aligned with one axis.
                            Array<int, 3> sumDirection[3];
                            sumDirection[0] =
                                inwardWallNormal[0] == 0
                                    ? Array<int, 3>(0, 0, 0)
                                    : (inwardWallNormal[0] > 0 ? Array<int, 3>(1, 0, 0)
                                                               : Array<int, 3>(-1, 0, 0));
                            sumDirection[1] =
                                inwardWallNormal[1] == 0
                                    ? Array<int, 3>(0, 0, 0)
                                    : (inwardWallNormal[1] > 0 ? Array<int, 3>(0, 1, 0)
                                                               : Array<int, 3>(0, -1, 0));
                            sumDirection[2] =
                                inwardWallNormal[2] == 0
                                    ? Array<int, 3>(0, 0, 0)
                                    : (inwardWallNormal[2] > 0 ? Array<int, 3>(0, 0, 1)
                                                               : Array<int, 3>(0, 0, -1));
                            T sum[3] = {
                                std::numeric_limits<T>::max(), std::numeric_limits<T>::max(),
                                std::numeric_limits<T>::max()};
                            for (int iSum = 0; iSum < 3; iSum++) {
                                if (sumDirection[iSum][0] + sumDirection[iSum][1]
                                        + sumDirection[iSum][2]
                                    != 0) {
                                    sum[iSum] = 0.0;
                                    for (int d = 0; d <= 3; d++) {
                                        plint posX = iX + d * sumDirection[iSum][0];
                                        plint posY = iY + d * sumDirection[iSum][1];
                                        plint posZ = iZ + d * sumDirection[iSum][2];
                                        if (!isAnyWall(param.flag(posX, posY, posZ))) {
                                            sum[iSum] += param.volumeFraction(posX, posY, posZ);
                                        }
                                    }
                                }
                            }

                            // The wall normal direction is the direction of the smallest sum.
                            if (sum[0] < sum[1]) {
                                if (sum[0] < sum[2]) {
                                    iWallNormalDirection = 0;
                                }
                                // sum[0]<sum[1] && sum[0] >= sum[2]
                                else
                                {
                                    iWallNormalDirection = 2;
                                }
                            }
                            // sum[0] >= sum[1]
                            else
                            {
                                if (sum[1] < sum[2]) {
                                    iWallNormalDirection = 1;
                                }
                                // sum[0] >= sum[1] && sum[1] >= sum[2]
                                else
                                {
                                    iWallNormalDirection = 2;
                                }
                            }
                        }

                        // Reset the inward wall normal to be unitary and to contain information on
                        // the direction.
                        inwardWallNormal[0] =
                            iWallNormalDirection != 0 ? 0 : (inwardWallNormal[0] > 0 ? 1 : -1);
                        inwardWallNormal[1] =
                            iWallNormalDirection != 1 ? 0 : (inwardWallNormal[1] > 0 ? 1 : -1);
                        inwardWallNormal[2] =
                            iWallNormalDirection != 2 ? 0 : (inwardWallNormal[2] > 0 ? 1 : -1);

                        // Define a wall normal that shows only the wall normal axis.
                        Array<int, 3> wallNormal;
                        wallNormal[0] = iWallNormalDirection == 0 ? 1 : 0;
                        wallNormal[1] = iWallNormalDirection == 1 ? 1 : 0;
                        wallNormal[2] = iWallNormalDirection == 2 ? 1 : 0;

                        // Compute the wall tangent vectors.
                        int iWallTangentDirection0 = iWallNormalDirection == 0     ? 1
                                                     : (iWallNormalDirection == 1) ? 2
                                                                                   : 0;
                        int iWallTangentDirection1 = iWallNormalDirection == 0     ? 2
                                                     : (iWallNormalDirection == 1) ? 0
                                                                                   : 1;
                        Array<int, 3> wallTangent0;
                        wallTangent0[0] = iWallTangentDirection0 == 0 ? 1 : 0;
                        wallTangent0[1] = iWallTangentDirection0 == 1 ? 1 : 0;
                        wallTangent0[2] = iWallTangentDirection0 == 2 ? 1 : 0;
                        Array<int, 3> wallTangent1;
                        wallTangent1[0] = iWallTangentDirection1 == 0 ? 1 : 0;
                        wallTangent1[1] = iWallTangentDirection1 == 1 ? 1 : 0;
                        wallTangent1[2] = iWallTangentDirection1 == 2 ? 1 : 0;

                        // Locally smooth the volume fraction to compute an estimate of the 2D
                        // normal.
                        T svfcp = param.smooth(
                            *param.volumeFractionP(), param.volumeFractionOffset(), iX, iY, iZ);
                        plint posX, posY, posZ;
                        posX = iX - wallTangent0[0];
                        posY = iY - wallTangent0[1];
                        posZ = iZ - wallTangent0[2];
                        T svf00 = !isAnyWall(param.flag(posX, posY, posZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), posX,
                                      posY, posZ)
                                                                           : svfcp;
                        posX = iX + wallTangent0[0];
                        posY = iY + wallTangent0[1];
                        posZ = iZ + wallTangent0[2];
                        T svf01 = !isAnyWall(param.flag(posX, posY, posZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), posX,
                                      posY, posZ)
                                                                           : svfcp;
                        posX = iX - wallTangent1[0];
                        posY = iY - wallTangent1[1];
                        posZ = iZ - wallTangent1[2];
                        T svf10 = !isAnyWall(param.flag(posX, posY, posZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), posX,
                                      posY, posZ)
                                                                           : svfcp;
                        posX = iX + wallTangent1[0];
                        posY = iY + wallTangent1[1];
                        posZ = iZ + wallTangent1[2];
                        T svf11 = !isAnyWall(param.flag(posX, posY, posZ)) ? param.smooth(
                                      *param.volumeFractionP(), param.volumeFractionOffset(), posX,
                                      posY, posZ)
                                                                           : svfcp;

                        // Compute a normalized 2D grad(VF) (inward-pointing 2D normal).
                        Array<T, 2> gradVF2D;
                        gradVF2D[0] = 0.5 * (svf01 - svf00);
                        gradVF2D[1] = 0.5 * (svf11 - svf10);
                        T norm_gradVF2D = norm(gradVF2D);
                        if (util::isZero(norm_gradVF2D)) {
                            param.curvature(iX, iY, iZ) = 0.0;
                            param.setNormal(iX, iY, iZ, zeroVector);
                            continue;
                        }
                        gradVF2D /= norm_gradVF2D;

                        T abs02D = std::fabs(gradVF2D[0]);
                        T abs12D = std::fabs(gradVF2D[1]);

                        int integrationDirection2D = 1;  // wallTangent1.
                        if (abs02D > abs12D) {
                            integrationDirection2D = 0;  // wallTangent0.
                        }

                        T h2D[3];
                        computeHeights2D(
                            param, wallTangent0, wallTangent1, integrationDirection2D, iX, iY, iZ,
                            h2D);

                        T dh2D = 0.5 * (h2D[2] - h2D[0]);

                        T sgn2D = -gradVF2D[integrationDirection2D] < 0.0 ? -1.0 : 1.0;
                        Array<T, 2> normal2D;
                        if (integrationDirection2D == 0)
                        {  // With respect to wallTangent0 and wallTangent1.
                            normal2D = Array<T, 2>(sgn2D, -dh2D);
                        } else {
                            normal2D = Array<T, 2>(-dh2D, sgn2D);
                        }
                        T norm_normal2D = norm(normal2D);
                        if (util::isZero(norm_normal2D)) {
                            param.curvature(iX, iY, iZ) = 0.0;
                            param.setNormal(iX, iY, iZ, zeroVector);
                            continue;
                        }

                        Array<T, 3> normal;  // 3D outward unit normal.
                        T wallNormalComponent = norm_normal2D / tanContactAngle;
                        normal[0] = normal2D[0] * wallTangent0[0] + normal2D[1] * wallTangent1[0]
                                    + wallNormalComponent * wallNormal[0];
                        normal[1] = normal2D[0] * wallTangent0[1] + normal2D[1] * wallTangent1[1]
                                    + wallNormalComponent * wallNormal[1];
                        normal[2] = normal2D[0] * wallTangent0[2] + normal2D[1] * wallTangent1[2]
                                    + wallNormalComponent * wallNormal[2];
                        T norm_normal = norm(normal);
                        if (util::isZero(norm_normal)) {
                            param.setNormal(iX, iY, iZ, zeroVector);
                        } else {
                            param.setNormal(iX, iY, iZ, normal / norm_normal);
                        }

                        // Now compute the curvature.
                        // First compute the 3D height functions.
                        int integrationDirection;
                        if (integrationDirection2D == 0) {
                            integrationDirection =
                                wallTangent0[0] != 0 ? 0 : (wallTangent0[1] != 0 ? 1 : 2);
                        } else {
                            integrationDirection =
                                wallTangent1[0] != 0 ? 0 : (wallTangent1[1] != 0 ? 1 : 2);
                        }
                        T h[3][3];
                        computeHeights3D(param, integrationDirection, iX, iY, iZ, h);

                        // Determine the orientation of the elements of h.
                        int iTangentDirection0 = integrationDirection == 0     ? 1
                                                 : (integrationDirection == 1) ? 2
                                                                               : 0;
                        int iTangentDirection1 = integrationDirection == 0     ? 2
                                                 : (integrationDirection == 1) ? 0
                                                                               : 1;
                        Array<int, 3> tangent0;
                        tangent0[0] = iTangentDirection0 == 0 ? 1 : 0;
                        tangent0[1] = iTangentDirection0 == 1 ? 1 : 0;
                        tangent0[2] = iTangentDirection0 == 2 ? 1 : 0;
                        Array<int, 3> tangent1;
                        tangent1[0] = iTangentDirection1 == 0 ? 1 : 0;
                        tangent1[1] = iTangentDirection1 == 1 ? 1 : 0;
                        tangent1[2] = iTangentDirection1 == 2 ? 1 : 0;

                        int i0 = -1;
                        int j0 = -1;
                        if (inwardWallNormal[0] == tangent0[0] && inwardWallNormal[1] == tangent0[1]
                            && inwardWallNormal[2] == tangent0[2])
                        {
                            i0 = 0;
                        } else if (
                            inwardWallNormal[0] == -tangent0[0]
                            && inwardWallNormal[1] == -tangent0[1]
                            && inwardWallNormal[2] == -tangent0[2])
                        {
                            i0 = 2;
                        } else if (
                            inwardWallNormal[0] == tangent1[0] && inwardWallNormal[1] == tangent1[1]
                            && inwardWallNormal[2] == tangent1[2])
                        {
                            j0 = 0;
                        } else if (
                            inwardWallNormal[0] == -tangent1[0]
                            && inwardWallNormal[1] == -tangent1[1]
                            && inwardWallNormal[2] == -tangent1[2])
                        {
                            j0 = 2;
                        } else {
                            PLB_ASSERT(false);
                        }

                        Array<T, 3> v1, v2;  // In the wallTangent0, wallTangent1 base.
                        v1[0] = std::fabs(normal2D[0]);
                        v1[1] = std::fabs(normal2D[1]);
                        v1[2] = 0.0;
                        if (integrationDirection2D == 0) {
                            v2[0] = 1.0;
                            v2[1] = 0.0;
                            v2[2] = 0.0;
                        } else {
                            v2[0] = 0.0;
                            v2[1] = 1.0;
                            v2[2] = 0.0;
                        }
                        T cosAlpha = std::cos(angleBetweenVectors(v1, v2));
                        T correction = 1.0 / (tanContactAngle * cosAlpha);

                        if (i0 != -1 && j0 == -1) {
                            for (int d = 0; d < 3; d++) {
                                h[i0][d] = h[1][d] + correction;
                            }
                        } else if (i0 == -1 && j0 != -1) {
                            for (int d = 0; d < 3; d++) {
                                h[d][j0] = h[d][1] + correction;
                            }
                        } else {
                            PLB_ASSERT(false);
                        }

                        T dh0 = 0.5 * (h[2][1] - h[0][1]);
                        T dh1 = 0.5 * (h[1][2] - h[1][0]);

                        T dh00 = h[2][1] - 2.0 * h[1][1] + h[0][1];
                        T dh11 = h[1][2] - 2.0 * h[1][1] + h[1][0];

                        T dh01 = 0.25 * (h[2][2] - h[2][0] - h[0][2] + h[0][0]);

                        T value = -(dh00 + dh11 + dh00 * dh1 * dh1 + dh11 * dh0 * dh0
                                    - 2.0 * dh01 * dh0 * dh1)
                                  / std::pow((T)1.0 + dh0 * dh0 + dh1 * dh1, (T)1.5);

                        param.curvature(iX, iY, iZ) = value;
                    }
                }
            }
        }

        delete interfaceFlag;
    }
}

/* *************** Class FreeSurfaceComputeCurvature3D ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceComputeCurvature3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    static T degToRad = (T)3.14159265358979323844L / (T)180;

    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    Dot3D location = param.absOffset();

    int smooth =
        0;  // 0 for no normal vector field smoothing.
            // 1 for smoothing the normal vector field once before imposing the contact angle.
            // 2 for smoothing the normal vector field once after imposing the contact angle.

    TensorField3D<T, 3> *normalForContactAngle = 0;
    Dot3D ofsNCA(0, 0, 0);
    plint w = 0;

    if (smooth == 0 || smooth == 2) {
        // Tensor field to hold a temporary vector field of unit normals. (Include also a w-cell
        // layer around "domain".)
        w = smooth == 0 ? 1 : 2;
        plint nx = domain.getNx() + 2 * w;
        plint ny = domain.getNy() + 2 * w;
        plint nz = domain.getNz() + 2 * w;
        normalForContactAngle = new TensorField3D<T, 3>(nx, ny, nz);
        ofsNCA = Dot3D(-domain.x0 + w, -domain.y0 + w, -domain.z0 + w);
        for (plint iX = domain.x0 - w; iX <= domain.x1 + w; ++iX) {
            plint i = iX + ofsNCA.x;
            for (plint iY = domain.y0 - w; iY <= domain.y1 + w; ++iY) {
                plint j = iY + ofsNCA.y;
                for (plint iZ = domain.z0 - w; iZ <= domain.z1 + w; ++iZ) {
                    plint k = iZ + ofsNCA.z;
                    normalForContactAngle->get(i, j, k) = param.getNormal(iX, iY, iZ);

                    param.curvature(iX, iY, iZ) = 0.0;
                }
            }
        }
    } else {
        // Tensor field to hold a temporary vector field of unit normals. (Include also a w-cell
        // layer around "domain".) This normal field is smoothed once (early tests show that if the
        // normal field is smoothed twice, then the results are pretty bad).
        w = 1;
        plint nx = domain.getNx() + 2 * w;
        plint ny = domain.getNy() + 2 * w;
        plint nz = domain.getNz() + 2 * w;
        normalForContactAngle = new TensorField3D<T, 3>(nx, ny, nz);
        ofsNCA = Dot3D(-domain.x0 + w, -domain.y0 + w, -domain.z0 + w);
        for (plint iX = domain.x0 - w; iX <= domain.x1 + w; ++iX) {
            plint i = iX + ofsNCA.x;
            for (plint iY = domain.y0 - w; iY <= domain.y1 + w; ++iY) {
                plint j = iY + ofsNCA.y;
                for (plint iZ = domain.z0 - w; iZ <= domain.z1 + w; ++iZ) {
                    plint k = iZ + ofsNCA.z;
                    normalForContactAngle->get(i, j, k) =
                        param.template lbmSmooth<T, descriptors::AdvectionDiffusionD3Q7Descriptor>(
                            *param.normalP(), param.normalOffset(), iX, iY, iZ);

                    param.curvature(iX, iY, iZ) = 0.0;
                }
            }
        }
    }

    // Enforce contact angles.
    if (useContactAngle) {
        for (plint iX = domain.x0 - w; iX <= domain.x1 + w; ++iX) {
            plint i = iX + ofsNCA.x;
            for (plint iY = domain.y0 - w; iY <= domain.y1 + w; ++iY) {
                plint j = iY + ofsNCA.y;
                for (plint iZ = domain.z0 - w; iZ <= domain.z1 + w; ++iZ) {
                    plint k = iZ + ofsNCA.z;

                    T curv = 0.0;

                    // If we want to operate on the "sharp free surface", which is defined by the
                    // lattice cells with the "interface" flag, then the next clause must be used.
                    // if (param.flag(iX, iY, iZ) == interface) {
                    // If we want a "diffuse interface view", then instead of the above clause, we
                    // need to have one of the following, which operate on a "thick free surface"
                    // that is an outcome of the smoothing proceedure (diffusion) performed when
                    // computing the field of unit normals. There are two ways to do this: one is to
                    // use the original normal ("getNormal"), which is already "thick", and the
                    // other way is to use the "normalForContactAngle" which is localy smoothed one
                    // more time if "smooth = 1".
                    if (!isAnyWall(param.flag(iX, iY, iZ))
                        && !util::isZero(normSqr(param.getNormal(iX, iY, iZ)))) {
                        // if (!isAnyWall(param.flag(iX, iY, iZ)) &&
                        // !util::isZero(normSqr(normalForContactAngle->get(i, j, k)))) {
                        int isaContactAngleCell = 0;
                        int numWallCells = 0;
                        // Computation of the inward-pointing wall normal (towards the fluid).
                        Array<int, 3> tmpWallNormal(0, 0, 0);
                        plint h = 1;  // This cannot be greater than 1 when w = 2!
                        for (plint dx = -h; dx <= h; dx++) {
                            for (plint dy = -h; dy <= h; dy++) {
                                for (plint dz = -h; dz <= h; dz++) {
                                    int flg = param.flag(iX + dx, iY + dy, iZ + dz);
                                    if (isWall(flg)) {  // Here we want only the no-slip walls, not
                                                        // the free-slip ones.
                                        tmpWallNormal += Array<int, 3>(-dx, -dy, -dz);
                                        numWallCells++;
                                    }
                                }
                            }
                        }
                        Array<T, 3> wallNormal;
                        if (numWallCells != 0) {
                            int norm2tmpWallNormal = tmpWallNormal[0] * tmpWallNormal[0]
                                                     + tmpWallNormal[1] * tmpWallNormal[1]
                                                     + tmpWallNormal[2] * tmpWallNormal[2];
                            if (norm2tmpWallNormal != 0) {
                                T tmpNormWallNormal = std::sqrt((T)norm2tmpWallNormal);
                                wallNormal[0] = (T)tmpWallNormal[0] / tmpNormWallNormal;
                                wallNormal[1] = (T)tmpWallNormal[1] / tmpNormWallNormal;
                                wallNormal[2] = (T)tmpWallNormal[2] / tmpNormWallNormal;
                                isaContactAngleCell = 1;
                            }
                        }
                        if (isaContactAngleCell) {
                            // Here again, there are two possible normals one can use: the original
                            // normal ("getNormal"), which is already "thick", and the
                            // "normalForContactAngle" which is localy smoothed one more time if
                            // "smooth = 1".
                            Array<T, 3> normal = param.getNormal(iX, iY, iZ);
                            // Array<T,3> normal = normalForContactAngle->get(i, j, k);
                            T cosPhi = dot(normal, wallNormal);
                            Array<T, 3> wallTangent = normal - cosPhi * wallNormal;
                            T normWallTangent = norm(wallTangent);
                            if (!util::isZero(normWallTangent)) {
                                wallTangent /= normWallTangent;
                            } else {
                                // This option to not enforce a contact angle when the interface
                                // normal is parallel to the wall normal must be revised in the
                                // framework of long thin films.
                                continue;
                            }

                            // Enforce the contact angle by taking at most two actions:
                            //
                            // 1) Add a penalty-like contribution to the curvature, to implicitly
                            // correct
                            //    towards the equilibrium contact angle. For the form of the
                            //    penalty-like term, see: Attar et al, "Lattice Boltzmann method for
                            //    dynamic wetting problems", Journal of Colloid and Interface
                            //    Science, 335 (2009) 84-93.
                            //
                            //    When operating in very low resolution, one should consider not
                            //    using this penalty term, for better accuracy.
                            //
                            // 2) Modify the free surface normal at the local interface cell, so
                            // that
                            //    the next time the curvature is computed, the equilibrium contact
                            //    angle is considered.
                            //
                            // Both of these actions impose the contact angle indirectly, through
                            // Young's law which relates the curvature, the surface tension and the
                            // pressure drop on the free surface.

                            T thetaEq =
                                (contactAngleFunction == 0
                                     ? contactAngle
                                     : contactAngleFunction(
                                           iX + location.x, iY + location.y, iZ + location.z)
                                           * degToRad);  // Equilibrium contact angle in radians.

                            // curv = cosPhi - std::cos(thetaEq); // At convergence, this term is
                            // zero.
                            normal =
                                std::cos(thetaEq) * wallNormal + std::sin(thetaEq) * wallTangent;

                            // This is yet another way, that does not try to impose the equilibrium
                            // contact angle, but an extrapolation of it, taking under consideration
                            // the current value of the dynamic contact angle (phi): theta = (T) 2 *
                            // thetaEq - phi. After convergence is achieved, then: phi = thetaEq,
                            // and theta = thetaEq.
                            //
                            // Another way is that one can impose the curvature as:
                            // curv = cosPhi - std::cos(thetaEq), but the normal as:
                            // normal = std::cos(theta) * wallNormal + std::sin(theta) *
                            // wallTangent, with: theta = (T) 2 * thetaEq - phi. This is a mixture
                            // of the two previously discussed methods. Sometimes it gives nice
                            // results (it enforces the contact angle more strongly).

                            // T phi = angleBetweenVectors(normal, wallNormal);
                            // T theta = (T) 2 * thetaEq - phi;

                            // curv = std::cos(theta) - std::cos(thetaEq); // At convergence, this
                            // term is zero. normal = std::cos(theta) * wallNormal + std::sin(theta)
                            // * wallTangent;

                            normalForContactAngle->get(i, j, k) = normal;
                        }
                    }

                    param.curvature(iX, iY, iZ) = curv;
                }
            }
        }
    }

    TensorField3D<T, 3> *normalForCurvature = 0;
    Dot3D ofsNC(0, 0, 0);

    if (smooth == 2) {
        // Tensor field to hold a temporary vector field of unit normals. (Include also a t-cell
        // layer around "domain".) This normal field is smoothed once (early tests show that if the
        // normal field is smoothed twice, then the results are pretty bad).
        plint t = 1;
        plint nx = domain.getNx() + 2 * t;
        plint ny = domain.getNy() + 2 * t;
        plint nz = domain.getNz() + 2 * t;
        normalForCurvature = new TensorField3D<T, 3>(nx, ny, nz);
        ofsNC = Dot3D(-domain.x0 + t, -domain.y0 + t, -domain.z0 + t);
        for (plint iX = domain.x0 - t; iX <= domain.x1 + t; ++iX) {
            plint i = iX + ofsNC.x;
            for (plint iY = domain.y0 - t; iY <= domain.y1 + t; ++iY) {
                plint j = iY + ofsNC.y;
                for (plint iZ = domain.z0 - t; iZ <= domain.z1 + t; ++iZ) {
                    plint k = iZ + ofsNC.z;
                    normalForCurvature->get(i, j, k) =
                        param.template lbmSmooth<T, descriptors::AdvectionDiffusionD3Q7Descriptor>(
                            *normalForContactAngle, ofsNCA, iX, iY, iZ);
                }
            }
        }
    } else {
        normalForCurvature = normalForContactAngle;
        ofsNC = ofsNCA;
    }

    // Compute the curvature as the divergence of the vector field of unit normals.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                T curv = 0.0;

                if (param.flag(iX, iY, iZ) != interface) {
                    param.curvature(iX, iY, iZ) = curv;
                    continue;
                }

                /*
                int useLB = 1;
                for (plint iPop = 1; iPop < Descriptor<T>::q; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    if (isAnyWall(param.flag(nextX, nextY, nextZ))) {
                        useLB = 0;
                        break;
                    }
                }

                if (useLB) {
                */
                // Compute the divergence of the normal vector field "the lattice Boltzmann way".
                curv = param.template lbmComputeDivergence<T, Descriptor>(
                    *normalForCurvature, ofsNC, iX, iY, iZ);
                /*
                } else {
                    // Compute the divergence of the normal vector field with finite differences on
                the interface cells
                    // excluding wall cells.
                    plint h = 1;    // This must be 1 because above we included only a 1-cell layer
                around "domain". curv = param.computeDivergence(*normalForCurvature, ofsNC, h, iX,
                iY, iZ);
                }
                */

                // We restrict the radius of curvature to be more always >= 0.5, in lattice units.
                // A smaller radius makes no sense anyway, numerically speaking, and in this way
                // we avoid problems of the "division by zero" kind. (radius = 2/curvature)
                // if (std::fabs(curv) > (T) 4) {
                //    if (curv < (T) 0) {
                //        curv = -4.0;
                //    }
                //    else {
                //        curv = 4.0;
                //    }
                //}

                param.curvature(iX, iY, iZ) +=
                    curv;  // We add in order to include the "penalty-like" term if any.
            }
        }
    }

    delete normalForContactAngle;
    if (smooth == 2) {
        delete normalForCurvature;
    }
}

/* *************** Class FreeSurfaceMassChange3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceMassChange3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    // This loop updates the mass, summarizing  Eq. 6/7, and Eq.8, in
    // the N. Thuerey e.a. technical report "Interactive Free Surface Fluids
    // with the Lattice Boltzmann Method".
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Cell<T, Descriptor> &cell = param.cell(iX, iY, iZ);
                int flag = param.flag(iX, iY, iZ);
                if (isFullWet(flag)) {
                    freeSurfaceTemplates<T, Descriptor>::massExchangeFluidCell(param, iX, iY, iZ);
                } else if (flag == interface) {
                    for (plint iPop = 1; iPop < D::q; ++iPop) {
                        plint nextX = iX + D::c[iPop][0];
                        plint nextY = iY + D::c[iPop][1];
                        plint nextZ = iZ + D::c[iPop][2];
                        int nextFlag = param.flag(nextX, nextY, nextZ);
                        plint opp = indexTemplates::opposite<D>(iPop);
                        // Calculate mass at time t+1 on interface cell --> eq 7 Thurey's paper.
                        if (isFullWet(nextFlag)) {
                            param.mass(iX, iY, iZ) +=
                                (cell[opp] - param.cell(nextX, nextY, nextZ)[iPop]);
                        } else if (nextFlag == interface) {
                            param.mass(iX, iY, iZ) +=
                                (cell[opp] - param.cell(nextX, nextY, nextZ)[iPop]) * 0.5
                                * (param.volumeFraction(nextX, nextY, nextZ)
                                   + param.volumeFraction(iX, iY, iZ));
                        }
                    }
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceCompletion3D ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceCompletion3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    using namespace freeSurfaceFlag;
    typedef typename InterfaceLists<T, Descriptor>::Node Node;

    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    // In this data processor, populations are both written locally and read non-locally.
    // To guarantee data consistency, a first loop makes only read accesses and stores
    // the necessary information into the list neighborOppositePop. A second loop reads
    // from this list and assigns values to populations.
    std::map<Node, Array<T, D::q> > neighborOppositePop;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                // This is the old form of the completion scheme. There is this extra condition
                // mentioned by Thurey which has to do with the normal to the interface. We found
                // that this condition is responsible for an instability when one increases
                // both the spatial and temporal resolution while respecting the diffusive limit
                // in the presence of surface tension. We also found that it causes an instability
                // at the simple test case of a fluid sphere which is subject to surface tension
                // but not to any other force. This sphere should remain still, but in the presence
                // of this condition it starts moving.
                /*
                if (param.flag(iX,iY,iZ) == interface) {
                    // Here we are on an interface node. The entire set of fi's is reconstructed.
                    // The normal is recomputed as in eq. 10 of Thurey's paper.
                    Array<T,3> normalToInterface;
                    normalToInterface = param.getNormal(iX, iY, iZ);

                    bool needsModification = false;
                    Array<T,D::q> savedPop;
                    savedPop[0] = -2.;
                    for(plint iPop=1; iPop < D::q; ++iPop )
                    {
                        // This is one of the tricky points of the code
                        // we have to decide if the f_is from the neighborhood
                        // have to be re-update by using the Thurey's rule, which
                        // states that f_i's coming from nearest neighs. that are empty cells,
                        // have to be re-updated.
                        // I like the eq.   f^{in}_i(x,t+dt) = f^{out}_i(x-e_i,t);
                        // This eq. makes me think that the neigh. that I have to check
                        // (to control is status e.g. empty or fluid ?) has to be pos-c_i
                        plint prevX = iX-D::c[iPop][0];
                        plint prevY = iY-D::c[iPop][1];
                        plint prevZ = iZ-D::c[iPop][2];

                        plint opp = indexTemplates::opposite<D>(iPop);
                        T scalarProduct = D::c[opp][0]*normalToInterface[0] +
                                          D::c[opp][1]*normalToInterface[1] +
                                          D::c[opp][2]*normalToInterface[2];

                        // Should I also change particle distribution function coming from
                        // bounceBack nodes? Well ideally no ... but there is for sure some
                        // cell configuration where these f_is are not well defined because
                        // they are probably coming from empty cells

                        // If the f_i[iPop] would be streamed from an empty cell, or whenever the
                scalar product is positive. if ( scalarProduct > 0 || param.flag(prevX,prevY,prevZ)
                == empty || isAnyWall(param.flag(prevX,prevY,prevZ)) )
                        {
                            savedPop[iPop] = param.cell(prevX,prevY,prevZ)[opp];
                            needsModification = true;
                        }
                        else {
                            savedPop[iPop] = (T)-2.;
                        }
                    }
                    if (needsModification) {
                        neighborOppositePop.insert(std::pair<Node,Array<T,D::q> >(Node(iX,iY,iZ),
                savedPop));
                    }
                }
                */

                if (param.flag(iX, iY, iZ) == interface) {
                    // Here we are on an interface node. The entire set of fi's is reconstructed.
                    bool needsModification = false;
                    Array<T, D::q> savedPop;
                    savedPop[0] = -2.;
                    for (plint iPop = 1; iPop < D::q; ++iPop) {
                        // This is one of the tricky points of the code
                        // we have to decide if the f_is from the neighborhood
                        // have to be re-update by using the Thurey's rule, which
                        // states that f_i's coming from nearest neighs. that are empty cells,
                        // have to be re-updated.
                        // I like the eq.   f^{in}_i(x,t+dt) = f^{out}_i(x-e_i,t);
                        // This eq. makes me think that the neigh. that I have to check
                        // (to control is status e.g. empty or fluid ?) has to be pos-c_i
                        plint prevX = iX - D::c[iPop][0];
                        plint prevY = iY - D::c[iPop][1];
                        plint prevZ = iZ - D::c[iPop][2];

                        plint opp = indexTemplates::opposite<D>(iPop);

                        // Should I also change particle distribution function coming from
                        // bounceBack nodes? Well ideally no ... but there is for sure some
                        // cell configuration where these f_is are not well defined because
                        // they are probably coming from empty cells

                        // If the f_i[iPop] would be streamed from an empty cell
                        if (isEmpty(param.flag(prevX, prevY, prevZ))
                            || isAnyWall(param.flag(prevX, prevY, prevZ))) {
                            savedPop[iPop] = param.cell(prevX, prevY, prevZ)[opp];
                            needsModification = true;
                        } else {
                            savedPop[iPop] = (T)-2.;
                        }
                    }
                    if (needsModification) {
                        neighborOppositePop.insert(
                            std::pair<Node, Array<T, D::q> >(Node(iX, iY, iZ), savedPop));
                    }
                }
            }
        }
    }

    typename std::map<Node, Array<T, D::q> >::const_iterator nodes = neighborOppositePop.begin();
    for (; nodes != neighborOppositePop.end(); ++nodes) {
        Node node = nodes->first;
        plint iX = node[0];
        plint iY = node[1];
        plint iZ = node[2];
        Array<T, D::q> neighborOppPop = nodes->second;
        for (plint iPop = 1; iPop < D::q; ++iPop) {
            if (neighborOppPop[iPop] > (T)-1.) {
                // Velocity is simply taken from the previous time step.
                Array<T, 3> j = param.getMomentum(iX, iY, iZ);
                T jSqr = VectorTemplate<T, Descriptor>::normSqr(j);
                // Remember: the value of pressure on an interface node has been set in
                // F, and is equal to the ambient pressure for a
                // single free-surface fluid, or in the case of a binary pressure, an
                // averaged value.
                T rhoBar = Descriptor<T>::rhoBar(param.getDensity(iX, iY, iZ));
                T feq_i = param.cell(iX, iY, iZ).computeEquilibrium(iPop, rhoBar, j, jSqr);
                plint opp = indexTemplates::opposite<D>(iPop);
                T feq_opp_i = param.cell(iX, iY, iZ).computeEquilibrium(opp, rhoBar, j, jSqr);
                param.cell(iX, iY, iZ)[iPop] = feq_i + feq_opp_i - neighborOppPop[iPop];
            }
        }
    }
}

/* *************** Class FreeSurfaceMacroscopic3D ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceMacroscopic3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    T lostMass = param.getSumLostMass();
    plint numInterfaceCells = param.getNumInterfaceCells();
    T massPerCell = T();
    if (numInterfaceCells > 0) {
        massPerCell = lostMass / (T)numInterfaceCells;
    }

    // Save macroscopic fields in external scalars and update the mass and the volume-fraction.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (isWet(param.flag(iX, iY, iZ))) {
                    T rhoBar;
                    Array<T, 3> j;
                    momentTemplates<T, Descriptor>::get_rhoBar_j(param.cell(iX, iY, iZ), rhoBar, j);
                    T density = Descriptor<T>::fullRho(rhoBar);
                    param.setDensity(iX, iY, iZ, density);

                    if (param.flag(iX, iY, iZ) == interface) {
                        param.mass(iX, iY, iZ) += massPerCell;
                        T newDensity = param.outsideDensity(iX, iY, iZ);
                        param.volumeFraction(iX, iY, iZ) = param.mass(iX, iY, iZ) / newDensity;
                        // On interface cells, adjust the pressure to the ambient pressure.
                        param.setDensity(iX, iY, iZ, newDensity);
                        if (!incompressibleModel) {
                            j *= newDensity / density;
                        }
                    } else if (isFullWet(param.flag(iX, iY, iZ))) {
                        param.volumeFraction(iX, iY, iZ) = T(1);
                    }

                    param.setMomentum(iX, iY, iZ, j);
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceMacroscopicWithoutLostMassReDistribution3D
 * ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceMacroscopicWithoutLostMassReDistribution3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    // Save macroscopic fields in external scalars and update the mass and the volume-fraction.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (isWet(param.flag(iX, iY, iZ))) {
                    T rhoBar;
                    Array<T, 3> j;
                    momentTemplates<T, Descriptor>::get_rhoBar_j(param.cell(iX, iY, iZ), rhoBar, j);
                    T density = Descriptor<T>::fullRho(rhoBar);
                    param.setDensity(iX, iY, iZ, density);

                    if (param.flag(iX, iY, iZ) == interface) {
                        T newDensity = param.outsideDensity(iX, iY, iZ);
                        param.volumeFraction(iX, iY, iZ) = param.mass(iX, iY, iZ) / newDensity;
                        // On interface cells, adjust the pressure to the ambient pressure.
                        param.setDensity(iX, iY, iZ, newDensity);
                        if (!incompressibleModel) {
                            j *= newDensity / density;
                        }
                    } else if (isFullWet(param.flag(iX, iY, iZ))) {
                        param.volumeFraction(iX, iY, iZ) = T(1);
                    }

                    param.setMomentum(iX, iY, iZ, j);
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceAddSurfaceTension3D ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceAddSurfaceTension3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    // Save macroscopic fields in external scalars and add the surface tension effect.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (param.flag(iX, iY, iZ) == interface) {
                    T density = param.getDensity(iX, iY, iZ);
                    T newDensity = density;
                    // Stored curvature is computed to be twice the mean curvature.
                    newDensity += surfaceTension * param.curvature(iX, iY, iZ) * D::invCs2;
                    param.volumeFraction(iX, iY, iZ) = param.mass(iX, iY, iZ) / newDensity;
                    // On interface cells, adjust the pressure to incorporate surface tension.
                    param.setDensity(iX, iY, iZ, newDensity);
                    if (!incompressibleModel) {
                        Array<T, 3> j = param.getMomentum(iX, iY, iZ);
                        Array<T, 3> newJ = j * newDensity / density;
                        param.setMomentum(iX, iY, iZ, newJ);
                    }
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceAddSurfaceTensionWeberModel3D ********************************
 */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceAddSurfaceTensionWeberModel3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    // Save macroscopic fields in external scalars and add the surface tension effect.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (param.flag(iX, iY, iZ) == interface) {
                    Cell<T, Descriptor> &cell = param.cell(iX, iY, iZ);
                    Array<T, 3> jTmp;
                    T rhoBarTmp;
                    cell.getDynamics().computeRhoBarJ(
                        cell, rhoBarTmp, jTmp);  // This does not include the free-surface BC.
                    T jTmpSqr = normSqr(jTmp);
                    T C = (T)0;
                    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                        T f = cell[iPop];
                        T feq =
                            cell.getDynamics().computeEquilibrium(iPop, rhoBarTmp, jTmp, jTmpSqr);
                        T fneq = f - feq;
                        C += std::fabs(fneq / fullF<T, Descriptor>(feq, iPop));
                    }
                    C /= (T)Descriptor<T>::q;

                    T density = param.getDensity(iX, iY, iZ);  // These include the free-surface BC.
                    Array<T, 3> j = param.getMomentum(iX, iY, iZ);
                    T newSurfaceTension = surfaceTension;
                    if (C <= criticalKn) {
                        T uSqr =
                            incompressibleModel ? normSqr(j) : normSqr(j) / (density * density);
                        T We = rhoDefault * uSqr * characteristicLength / surfaceTension;
                        if (We > criticalWe) {
                            newSurfaceTension *= We / criticalWe;
                        }
                    }

                    T newDensity = density;
                    // Stored curvature is computed to be twice the mean curvature.
                    newDensity += newSurfaceTension * param.curvature(iX, iY, iZ) * D::invCs2;
                    param.volumeFraction(iX, iY, iZ) = param.mass(iX, iY, iZ) / newDensity;
                    // On interface cells, adjust the pressure to incorporate surface tension.
                    param.setDensity(iX, iY, iZ, newDensity);
                    if (!incompressibleModel) {
                        Array<T, 3> newJ = j * newDensity / density;
                        param.setMomentum(iX, iY, iZ, newJ);
                    }
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceAddSurfaceTensionFromScalarField3D
 * ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceAddSurfaceTensionFromScalarField3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    using namespace freeSurfaceFlag;

    std::vector<AtomicBlock3D *> freeSurfaceBlocks(atomicBlocks.begin(), atomicBlocks.end() - 1);
    FreeSurfaceProcessorParam3D<T, Descriptor> param(freeSurfaceBlocks);

    ScalarField3D<T> *surfaceTensionField =
        dynamic_cast<ScalarField3D<T> *>(*(atomicBlocks.end() - 1));
    PLB_ASSERT(surfaceTensionField);

    Dot3D ofsSTF = computeRelativeDisplacement(*atomicBlocks[0], *surfaceTensionField);

    // Save macroscopic fields in external scalars and add the surface tension effect.
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (param.flag(iX, iY, iZ) == interface) {
                    T density = param.getDensity(iX, iY, iZ);
                    T newDensity = density;
                    T surfaceTension =
                        surfaceTensionField->get(iX + ofsSTF.x, iY + ofsSTF.y, iZ + ofsSTF.z);
                    // Stored curvature is computed to be twice the mean curvature.
                    newDensity += surfaceTension * param.curvature(iX, iY, iZ) * D::invCs2;
                    param.volumeFraction(iX, iY, iZ) = param.mass(iX, iY, iZ) / newDensity;
                    // On interface cells, adjust the pressure to incorporate surface tension.
                    param.setDensity(iX, iY, iZ, newDensity);
                    if (!incompressibleModel) {
                        Array<T, 3> j = param.getMomentum(iX, iY, iZ);
                        Array<T, 3> newJ = j * newDensity / density;
                        param.setMomentum(iX, iY, iZ, newJ);
                    }
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceStabilize3D ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceStabilize3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (param.flag(iX, iY, iZ) == interface) {
                    Cell<T, Descriptor> &cell = param.cell(iX, iY, iZ);
                    T oldRhoBar;
                    Array<T, 3> oldJ;
                    momentTemplates<T, Descriptor>::get_rhoBar_j(cell, oldRhoBar, oldJ);
                    T oldJsqr = normSqr(oldJ);
                    T newDensity = param.getDensity(iX, iY, iZ);
                    T newRhoBar = Descriptor<T>::rhoBar(newDensity);
                    Array<T, 3> newJ = param.getMomentum(iX, iY, iZ);
                    T newJsqr = normSqr(newJ);
                    for (int iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                        T oldEq =
                            cell.getDynamics().computeEquilibrium(iPop, oldRhoBar, oldJ, oldJsqr);
                        T newEq =
                            cell.getDynamics().computeEquilibrium(iPop, newRhoBar, newJ, newJsqr);
                        cell[iPop] += newEq - oldEq;
                    }
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceComputeInterfaceLists3D
 * ******************************************* */

template <typename T, template <typename> class Descriptor>
T FreeSurfaceComputeInterfaceLists3D<T, Descriptor>::kappa = 1.e-3;

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceComputeInterfaceLists3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T, Descriptor>::Node Node;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);
    using namespace freeSurfaceFlag;

    param.massExcess().clear();
    param.interfaceToFluid().clear();
    param.interfaceToEmpty().clear();
    param.emptyToInterface().clear();

    // interfaceToFluid needs to be computed in bulk+2.
    for (plint iX = domain.x0 - 2; iX <= domain.x1 + 2; ++iX) {
        for (plint iY = domain.y0 - 2; iY <= domain.y1 + 2; ++iY) {
            for (plint iZ = domain.z0 - 2; iZ <= domain.z1 + 2; ++iZ) {
                Node node(iX, iY, iZ);
                // Eq. 11 in Thuerey's technical report.
                if (param.flag(iX, iY, iZ) == interface) {  // Interface cell.
                    if (param.volumeFraction(iX, iY, iZ) > T(1) + kappa)
                    {  // Interface cell is filled.
                        // Elements are added even if they belong to the envelope, because they may
                        // be
                        //   needed further down in the same data processor.
                        param.interfaceToFluid().insert(node);
                    } else if (param.volumeFraction(iX, iY, iZ) < kappa)
                    {  // Interface cell is empty.
                        // Elements are added even if they belong to the envelope, because they may
                        // be
                        //   needed further down in the same data processor.
                        param.interfaceToEmpty().insert(node);
                    }
                }
            }
        }
    }

    // Where interface cells have become fluid, neighboring cells must be prevented from
    //   being empty, because otherwise there's no interface cell between empty and fluid.
    typename std::set<Node>::iterator iEle = param.interfaceToFluid().begin();
    for (; iEle != param.interfaceToFluid().end(); ++iEle) {
        // The node here may belong to the 1st envelope.
        Node node = *iEle;
        plint iX = node[0];
        plint iY = node[1];
        plint iZ = node[2];

        for (plint iPop = 1; iPop < D::q; ++iPop) {
            plint nextX = iX + D::c[iPop][0];
            plint nextY = iY + D::c[iPop][1];
            plint nextZ = iZ + D::c[iPop][2];
            Node nextNode(nextX, nextY, nextZ);

            // If one of my neighbors switches interface->fluid, then I shall be prevented
            //     from switching interface->empty at the same time step.
            if (contained(nextX, nextY, nextZ, domain.enlarge(1))
                && param.flag(nextX, nextY, nextZ) == interface)
            {
                param.interfaceToEmpty().erase(nextNode);
            }
            // If one of my neighbors switches interface->fluid and I am empty I shall become
            //   interface.
            else if (
                contained(nextX, nextY, nextZ, domain.enlarge(1))
                && isEmpty(param.flag(nextX, nextY, nextZ)))
            {
                param.emptyToInterface().insert(nextNode);
            }
        }
    }
}

/* *************** Class FreeSurfaceIniInterfaceToAnyNodes3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceIniInterfaceToAnyNodes3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T, Descriptor>::Node Node;
    using namespace freeSurfaceFlag;

    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    // 1. For interface->fluid nodes, update in the flag matrix,
    //   and compute and store mass excess from these cells.
    typename std::set<Node>::iterator iEle = param.interfaceToFluid().begin();
    for (; iEle != param.interfaceToFluid().end(); ++iEle) {
        Node node = *iEle;

        plint iX = node[0];
        plint iY = node[1];
        plint iZ = node[2];

        if (contained(iX, iY, iZ, domain.enlarge(1))) {
            T saveMass = param.mass(iX, iY, iZ);
            param.mass(iX, iY, iZ) = param.getDensity(iX, iY, iZ);
            param.volumeFraction(iX, iY, iZ) = (T)1;
            param.flag(iX, iY, iZ) = fluid;

            T massExcess = saveMass - param.getDensity(iX, iY, iZ);
            param.massExcess().insert(std::pair<Node, T>(node, massExcess));
        }
    }

    // 2. For interface->empty nodes, update in the flag matrix,
    //   and compute and store mass excess from these cells.
    iEle = param.interfaceToEmpty().begin();
    for (; iEle != param.interfaceToEmpty().end(); ++iEle) {
        Node node = *iEle;
        plint iX = node[0];
        plint iY = node[1];
        plint iZ = node[2];

        if (contained(iX, iY, iZ, domain.enlarge(1))) {
            // Avoid the case where an empty cell has a fluid neighbor without
            // interface cell between them.
            bool isAdjacentToProtected = false;
            for (plint iPop = 1; iPop < D::q; ++iPop) {
                plint nextX = iX + D::c[iPop][0];
                plint nextY = iY + D::c[iPop][1];
                plint nextZ = iZ + D::c[iPop][2];
                if (isProtected(param.flag(nextX, nextY, nextZ))) {
                    isAdjacentToProtected = true;
                    break;
                }
            }
            if (!isAdjacentToProtected) {
                param.flag(iX, iY, iZ) = empty;
                param.attributeDynamics(iX, iY, iZ, emptyNodeDynamicsTemplate->clone());

                T massExcess = param.mass(iX, iY, iZ);
                param.massExcess().insert(std::pair<Node, T>(node, massExcess));

                param.mass(iX, iY, iZ) = T();
                param.volumeFraction(iX, iY, iZ) = T();
                param.setDensity(iX, iY, iZ, rhoDefault);
                param.setForce(iX, iY, iZ, Array<T, D::ExternalField::sizeOfForce>::zero());
                param.setMomentum(iX, iY, iZ, Array<T, 3>(T(), T(), T()));
                for (plint iPop = 1; iPop < D::q; ++iPop) {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    plint nextZ = iZ + D::c[iPop][2];

                    // The concurrent read/write on param.flag is not an issue here, because the
                    // result in any case is that all adjacent fluid cells have become interface.
                    if (param.flag(nextX, nextY, nextZ) == fluid) {
                        param.flag(nextX, nextY, nextZ) = interface;
                    }
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceIniEmptyToInterfaceNodes3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceIniEmptyToInterfaceNodes3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T, Descriptor>::Node Node;
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    // In this data processor, density and momentum are potentially read and written
    //   from the same node, because nodes can switch state. The following two vectors
    //   store temporary variables to avoid read/write in undefined order.
    std::vector<T> newDensity(param.emptyToInterface().size());
    std::vector<Array<T, 3> > newMomentum(param.emptyToInterface().size());
    std::fill(newDensity.begin(), newDensity.end(), T());
    std::fill(newMomentum.begin(), newMomentum.end(), Array<T, 3>(T(), T(), T()));

    // Compute density and momentum for cells that will switch state empty->interface.
    //   It is sufficient to do this is bulk+0.
    //   This loop performs read-only access to the lattice.
    plint i = 0;
    typename std::set<Node>::iterator iEle = param.emptyToInterface().begin();
    for (; iEle != param.emptyToInterface().end(); ++iEle, ++i) {
        Node node = *iEle;
        plint iX = node[0];
        plint iY = node[1];
        plint iZ = node[2];

        // If non-bulk elements are left in the list, disregard to avoid accessing undefined
        // neighbors.
        if (contained(iX, iY, iZ, domain)) {
            // For initialization of the new cell, compute average density
            //   and momentum on neighbors.
            T averageDensity = T(0);
            Array<T, 3> averageMomentum(T(0), T(0), T(0));
            T sumWeights = (T)0;
            for (plint iPop = 1; iPop < D::q; ++iPop) {
                plint nextX = iX + D::c[iPop][0];
                plint nextY = iY + D::c[iPop][1];
                plint nextZ = iZ + D::c[iPop][2];

                // Warning: it is not accounted for the fact that neighbors can have excess mass. It
                //   might be good to account for this in the future.
                if (isWet(param.flag(nextX, nextY, nextZ))) {
                    T weight = D::t[iPop];
                    sumWeights += weight;
                    averageDensity += weight * param.getDensity(nextX, nextY, nextZ);
                    averageMomentum += weight * param.getMomentum(nextX, nextY, nextZ);
                }
            }
            T invSum = T(1) / sumWeights;
            averageDensity *= invSum;
            averageMomentum *= invSum;
            newDensity[i] = averageDensity;
            newMomentum[i] = averageMomentum;
        }
    }

    // Elements that have switched state empty->interface are initialized at equilibrium.
    //   It is sufficient to initialize them in bulk+0.
    //   This loop performs write-only access on the lattice.
    i = 0;
    iEle = param.emptyToInterface().begin();
    for (; iEle != param.emptyToInterface().end(); ++iEle, ++i) {
        Node node = *iEle;

        plint iX = node[0];
        plint iY = node[1];
        plint iZ = node[2];

        // If non-bulk elements are left in the list, disregard to avoid accessing undefined
        // neighbors.
        if (contained(iX, iY, iZ, domain)) {
            T averageDensity = newDensity[i];
            Array<T, 3> averageMomentum = newMomentum[i];

            param.attributeDynamics(iX, iY, iZ, dynamicsTemplate->clone());

            iniCellAtEquilibrium(
                param.cell(iX, iY, iZ), averageDensity, averageMomentum / averageDensity);
            param.setForce(iX, iY, iZ, force);
            // Change density, but leave mass and volumeFraction at 0, as they are later
            //   recomputed (Warning: this is probably correct, but there remains a small doubt).
            param.setMomentum(iX, iY, iZ, averageMomentum);
            param.setDensity(iX, iY, iZ, averageDensity);
            param.mass(iX, iY, iZ) = T();
            param.volumeFraction(iX, iY, iZ) = T();
            param.flag(iX, iY, iZ) = interface;
        }
    }
}

/* *************** Class FreeSurfaceRemoveFalseInterfaceCells3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceRemoveFalseInterfaceCells3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T, Descriptor>::Node Node;
    using namespace freeSurfaceFlag;

    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    /// In the following, the flag status of cells is read (non-locally) and
    /// modified (locally). To avoid conflict, two loops are made, the first
    /// of which reads only, and the second writes. The vectors "interfaceToFluidNodes"
    /// and "interfaceToEmptyNodes" store coordinates of nodes that will switch
    /// status.
    std::vector<Node> interfaceToFluidNodes, interfaceToEmptyNodes;
    for (plint iX = domain.x0 - 1; iX <= domain.x1 + 1; ++iX) {
        for (plint iY = domain.y0 - 1; iY <= domain.y1 + 1; ++iY) {
            for (plint iZ = domain.z0 - 1; iZ <= domain.z1 + 1; ++iZ) {
                Node node(iX, iY, iZ);
                if (param.flag(iX, iY, iZ) == interface) {
                    bool noEmptyNeighbor = true;
                    bool noFluidNeighbor = true;

                    for (plint iPop = 1; iPop < D::q; iPop++) {
                        plint nextX = iX + D::c[iPop][0];
                        plint nextY = iY + D::c[iPop][1];
                        plint nextZ = iZ + D::c[iPop][2];

                        if (isEmpty(param.flag(nextX, nextY, nextZ)))
                            noEmptyNeighbor = false;
                        if (isFullWet(param.flag(nextX, nextY, nextZ)))
                            noFluidNeighbor = false;
                    }

                    if (noEmptyNeighbor) {
                        interfaceToFluidNodes.push_back(Node(iX, iY, iZ));
                        // Store the coordinates, so flag on this node
                        // can be changed in a loop outside the current one.

                        T massExcess = param.mass(iX, iY, iZ) - param.getDensity(iX, iY, iZ);
                        param.massExcess().insert(std::pair<Node, T>(node, massExcess));
                        param.mass(iX, iY, iZ) = param.getDensity(iX, iY, iZ);
                        param.volumeFraction(iX, iY, iZ) = T(1);
                    } else if (noFluidNeighbor) {
                        interfaceToEmptyNodes.push_back(Node(iX, iY, iZ));
                        // Store the coordinates, so flag on this node
                        // can be changed in a loop outside the current one.

                        T massExcess = param.mass(iX, iY, iZ);
                        param.massExcess().insert(std::pair<Node, T>(node, massExcess));

                        param.attributeDynamics(iX, iY, iZ, emptyNodeDynamicsTemplate->clone());
                        param.mass(iX, iY, iZ) = T();
                        param.setDensity(iX, iY, iZ, rhoDefault);
                        param.volumeFraction(iX, iY, iZ) = T();
                        param.setForce(iX, iY, iZ, Array<T, D::ExternalField::sizeOfForce>::zero());
                        param.setMomentum(iX, iY, iZ, Array<T, 3>(T(), T(), T()));
                    }
                }
            }
        }
    }

    for (pluint i = 0; i < interfaceToFluidNodes.size(); ++i) {
        Node const &pos = interfaceToFluidNodes[i];
        param.flag(pos[0], pos[1], pos[2]) = fluid;
    }
    for (pluint i = 0; i < interfaceToEmptyNodes.size(); ++i) {
        Node const &pos = interfaceToEmptyNodes[i];
        param.flag(pos[0], pos[1], pos[2]) = empty;
    }
}

/* *************** Class FreeSurfaceEqualMassExcessReDistribution3D
 * ******************************************* */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceEqualMassExcessReDistribution3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T, Descriptor>::Node Node;
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    Box3D originalDomain(domain);

    typename std::map<Node, T>::iterator iEle = param.massExcess().begin();
    for (; iEle != param.massExcess().end(); ++iEle) {
        Array<plint, 3> node = iEle->first;
        plint iX = node[0];
        plint iY = node[1];
        plint iZ = node[2];

        // Check for valid interface neighbors to re-distribute mass
        if (contained(iX, iY, iZ, domain.enlarge(1))) {
            std::vector<int> indX, indY, indZ;
            plint numValidNeighbors = 0;

            // Check for interface neighbors in the LB directions.
            for (plint iPop = 1; iPop < D::q; iPop++) {
                plint nextX = iX + D::c[iPop][0];
                plint nextY = iY + D::c[iPop][1];
                plint nextZ = iZ + D::c[iPop][2];
                if (param.flag(nextX, nextY, nextZ) == interface) {
                    if (contained(nextX, nextY, nextZ, domain)) {
                        indX.push_back(nextX);
                        indY.push_back(nextY);
                        indZ.push_back(nextZ);
                    }
                    numValidNeighbors++;
                }
            }

            // Mass re-distribution
            if (numValidNeighbors != 0) {
                int indSize = (int)indX.size();
                T massToRedistribute = iEle->second / (T)numValidNeighbors;

                for (int i = 0; i < indSize; i++) {
                    int nextX = indX[i];
                    int nextY = indY[i];
                    int nextZ = indZ[i];

                    param.mass(nextX, nextY, nextZ) += massToRedistribute;
                    param.volumeFraction(nextX, nextY, nextZ) =
                        param.mass(nextX, nextY, nextZ) / param.getDensity(nextX, nextY, nextZ);
                }
            } else {
                if (contained(iX, iY, iZ, originalDomain)) {
                    param.addToLostMass(iEle->second);
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceEqualMassExcessReDistributionAndComputationOfLostMass3D
 * ********************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceEqualMassExcessReDistributionAndComputationOfLostMass3D<
    T, Descriptor>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    typedef Descriptor<T> D;
    typedef typename InterfaceLists<T, Descriptor>::Node Node;
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    Box3D originalDomain(domain);

    typename std::map<Node, T>::iterator iEle = param.massExcess().begin();
    for (; iEle != param.massExcess().end(); ++iEle) {
        Array<plint, 3> node = iEle->first;
        plint iX = node[0];
        plint iY = node[1];
        plint iZ = node[2];

        // Check for valid interface neighbors to re-distribute mass
        if (contained(iX, iY, iZ, domain.enlarge(1))) {
            std::vector<int> indX, indY, indZ;
            plint numValidNeighbors = 0;

            // Check for interface neighbors in the LB directions.
            for (plint iPop = 1; iPop < D::q; iPop++) {
                plint nextX = iX + D::c[iPop][0];
                plint nextY = iY + D::c[iPop][1];
                plint nextZ = iZ + D::c[iPop][2];
                if (param.flag(nextX, nextY, nextZ) == interface) {
                    if (contained(nextX, nextY, nextZ, domain)) {
                        indX.push_back(nextX);
                        indY.push_back(nextY);
                        indZ.push_back(nextZ);
                    }
                    numValidNeighbors++;
                }
            }

            // Mass re-distribution
            if (numValidNeighbors != 0) {
                int indSize = (int)indX.size();
                T massToRedistribute = iEle->second / (T)numValidNeighbors;

                for (int i = 0; i < indSize; i++) {
                    int nextX = indX[i];
                    int nextY = indY[i];
                    int nextZ = indZ[i];

                    param.mass(nextX, nextY, nextZ) += massToRedistribute;
                    param.volumeFraction(nextX, nextY, nextZ) =
                        param.mass(nextX, nextY, nextZ) / param.getDensity(nextX, nextY, nextZ);
                }
            } else {
                if (contained(iX, iY, iZ, originalDomain)) {
                    reductionData.localLostMass += iEle->second;
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceComputeStatistics3D *******************************************
 */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceComputeStatistics3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;

    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (isWet(param.flag(iX, iY, iZ))) {
                    param.addToTotalMass(param.mass(iX, iY, iZ));
                    if (param.flag(iX, iY, iZ) == interface) {
                        param.addToInterfaceCells(1);
                    }
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceComputeReductionsPerProcess3D
 * ******************************************* */

template <typename T>
void FreeSurfaceComputeReductionsPerProcess3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    using namespace freeSurfaceFlag;

    PLB_ASSERT(blocks.size() == 2);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int> *>(blocks[0]);
    PLB_ASSERT(flag);
    ScalarField3D<T> *mass = dynamic_cast<ScalarField3D<T> *>(blocks[1]);
    PLB_ASSERT(mass);

    Dot3D ofs = computeRelativeDisplacement(*flag, *mass);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (isWet(flag->get(iX, iY, iZ))) {
                    reductionData.localTotalMass += mass->get(iX + ofs.x, iY + ofs.y, iZ + ofs.z);
                    if (flag->get(iX, iY, iZ) == interface) {
                        reductionData.localNumInterfaceCells += 1;
                    }
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceAddExternalForce3D ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceAddExternalForce3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    if (Descriptor<T>::ExternalField::sizeOfForce == 0) {
        return;
    }
    PLB_ASSERT(Descriptor<T>::ExternalField::sizeOfForce == 3);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (isWet(param.flag(iX, iY, iZ))) {
                    Array<T, 3> newJ = param.getMomentum(iX, iY, iZ);
                    Array<T, Descriptor<T>::ExternalField::sizeOfForce> force =
                        param.getForce(iX, iY, iZ);

                    T dynamicOmega =
                        param.cell(iX, iY, iZ)
                            .getDynamics()
                            .getDynamicParameter(
                                dynamicParams::dynamicOmega,
                                param.cell(iX, iY, iZ));  // In case of a Smagorinsky model.
                    T tau = 0.0;
                    if (!util::isZero(dynamicOmega)) {
                        tau = (T)1 / dynamicOmega;
                    } else {
                        tau = (T)1 / param.cell(iX, iY, iZ).getDynamics().getOmega();
                    }

                    // Two comments:
                    // - Here the force is multiplied by rho0 and not rho so that, under
                    //   gravity, a linear pressure profile is obtained.
                    // - The force is not multiplied by the volume fraction (some authors
                    //   do multiply it by the volumeFraction), because there is a
                    //   point-wise interpretation of quantities like momentum.
                    for (plint i = 0; i < Descriptor<T>::ExternalField::sizeOfForce; i++) {
                        newJ[i] += rhoDefault * tau * force[i];
                    }
                    param.setMomentum(iX, iY, iZ, newJ);
                }
            }
        }
    }
}

/* *************** Class FreeSurfaceLostMassReDistribution3D ******************************** */

template <typename T, template <typename U> class Descriptor>
void FreeSurfaceLostMassReDistribution3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;
    FreeSurfaceProcessorParam3D<T, Descriptor> param(atomicBlocks);

    T lostMass = reductionData.lostMass;
    plint numInterfaceCells = reductionData.numInterfaceCells;
    T massPerCell = T();
    if (numInterfaceCells > 0) {
        massPerCell = lostMass / (T)numInterfaceCells;
    }

    if (massPerCell != (T)0) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    if (param.flag(iX, iY, iZ) == interface) {
                        param.mass(iX, iY, iZ) += massPerCell;
                    }
                }
            }
        }
    }
}

/* *************** Class RepelInterfaceFromImmersedWalls3D ******************************** */

template <typename T, class VelFunction>
void RepelInterfaceFromImmersedWalls3D<T, VelFunction>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;

    PLB_PRECONDITION(atomicBlocks.size() == 4);
    ScalarField3D<T> *rhoBar = dynamic_cast<ScalarField3D<T> *>(atomicBlocks[0]);
    TensorField3D<T, 3> *j = dynamic_cast<TensorField3D<T, 3> *>(atomicBlocks[1]);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int> *>(atomicBlocks[2]);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[3]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(j);
    PLB_ASSERT(flag);
    PLB_ASSERT(container);
    Dot3D location = rhoBar->getLocation();
    Dot3D ofsJ = computeRelativeDisplacement(*rhoBar, *j);
    Dot3D ofsF = computeRelativeDisplacement(*rhoBar, *flag);
    ImmersedWallData3D<T> *wallData = dynamic_cast<ImmersedWallData3D<T> *>(container->getData());
    PLB_ASSERT(wallData);

    std::vector<Array<T, 3> > const &vertices = wallData->vertices;
    std::vector<pluint> const &globalVertexIds = wallData->globalVertexIds;
    PLB_ASSERT(vertices.size() == globalVertexIds.size());

    if (strongRepelling) {
        // Define temporary matrices to facilitate the computation of the
        // mean position and velocity of the immersed boundary vertices
        // in a vicinity of each interface node.
        plint nx = domain.getNx();
        plint ny = domain.getNy();
        plint nz = domain.getNz();
        TensorField3D<T, 3> sumVelocity(nx, ny, nz);
        ScalarField3D<plint> numVertices(nx, ny, nz);
        Dot3D ofs(-domain.x0, -domain.y0, -domain.z0);

        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 3> const &vertex = vertices[i];
            Array<plint, 3> intPos(
                (plint)vertex[0] - location.x, (plint)vertex[1] - location.y,
                (plint)vertex[2] - location.z);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> zLim(
                (vertex[2] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
                for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                    for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                        Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                        if (contained(pos[0], pos[1], pos[2], domain)) {
                            if (flag->get(pos[0] + ofsF.x, pos[1] + ofsF.y, pos[2] + ofsF.z)
                                == interface) {
                                sumVelocity.get(pos[0] + ofs.x, pos[1] + ofs.y, pos[2] + ofs.z) +=
                                    velFunction(globalVertexIds[i]);
                                numVertices.get(pos[0] + ofs.x, pos[1] + ofs.y, pos[2] + ofs.z) +=
                                    1;
                            }
                        }
                    }
                }
            }
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    if (flag->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z)
                            == interface && numVertices.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z) != 0)
                    {
                        Array<T, 3> meanVelocityWall =
                            sumVelocity.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z)
                            / (T)numVertices.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z);

                        rhoBar->get(iX, iY, iZ) = rhoDefault - (T)1.0;
                        j->get(iX + ofsJ.x, iY + ofsJ.y, iZ + ofsJ.z) =
                            rhoDefault * meanVelocityWall;
                    }
                }
            }
        }
    } else {
        // Define temporary matrices to facilitate the computation of the
        // mean position and velocity of the immersed boundary vertices
        // in a vicinity of each interface node.
        plint nx = domain.getNx();
        plint ny = domain.getNy();
        plint nz = domain.getNz();
        TensorField3D<T, 3> sumPosition(nx, ny, nz);
        TensorField3D<T, 3> sumVelocity(nx, ny, nz);
        ScalarField3D<plint> numVertices(nx, ny, nz);
        Dot3D ofs(-domain.x0, -domain.y0, -domain.z0);

        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 3> const &vertex = vertices[i];
            Array<plint, 3> intPos(
                (plint)vertex[0] - location.x, (plint)vertex[1] - location.y,
                (plint)vertex[2] - location.z);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> zLim(
                (vertex[2] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
                for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                    for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                        Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                        if (contained(pos[0], pos[1], pos[2], domain)) {
                            if (flag->get(pos[0] + ofsF.x, pos[1] + ofsF.y, pos[2] + ofsF.z)
                                == interface) {
                                sumPosition.get(pos[0] + ofs.x, pos[1] + ofs.y, pos[2] + ofs.z) +=
                                    vertex;
                                sumVelocity.get(pos[0] + ofs.x, pos[1] + ofs.y, pos[2] + ofs.z) +=
                                    velFunction(globalVertexIds[i]);
                                numVertices.get(pos[0] + ofs.x, pos[1] + ofs.y, pos[2] + ofs.z) +=
                                    1;
                            }
                        }
                    }
                }
            }
        }

        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    if (flag->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z)
                            == interface && numVertices.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z) != 0)
                    {
                        // For each interface node which is close to an immersed boundary,
                        // compute the mean position vector of its neighbor empty nodes.
                        // This way we can calculate on which side of the immersed wall
                        // the specific bubble (or droplet) is situated.
                        Array<T, 3> meanPositionEmpty((T)0, (T)0, (T)0);
                        plint numEmpty = 0;
                        for (plint dx = -2; dx <= +2; ++dx) {
                            plint x = iX + dx;
                            for (plint dy = -2; dy <= +2; ++dy) {
                                plint y = iY + dy;
                                for (plint dz = -2; dz <= +2; ++dz) {
                                    plint z = iZ + dz;
                                    if (isEmpty(flag->get(x + ofsF.x, y + ofsF.y, z + ofsF.z))) {
                                        meanPositionEmpty += Array<T, 3>(
                                            (T)(x + location.x), (T)(y + location.y),
                                            (T)(z + location.z));
                                        numEmpty++;
                                    }
                                }
                            }
                        }
                        PLB_ASSERT(numEmpty);
                        meanPositionEmpty /= (T)numEmpty;

                        // Decide if the immersed boundary and the bubble (or droplet) are
                        // approaching one another. In such a case change the density and the
                        // momentum of the bubble (or droplet) interface nodes.
                        Array<T, 3> meanPositionWall =
                            sumPosition.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z)
                            / (T)numVertices.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z);
                        Array<T, 3> meanVelocityWall =
                            sumVelocity.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z)
                            / (T)numVertices.get(iX + ofs.x, iY + ofs.y, iZ + ofs.z);
                        T dotProd = dot(meanVelocityWall, meanPositionEmpty - meanPositionWall);

                        if (util::greaterEqual(dotProd, (T)0)) {
                            rhoBar->get(iX, iY, iZ) = rhoDefault - (T)1.0;
                            j->get(iX + ofsJ.x, iY + ofsJ.y, iZ + ofsJ.z) =
                                rhoDefault * meanVelocityWall;
                        }
                    }
                }
            }
        }
    }
}

/* *************** Class TemporarilyProtectImmersedWalls3D ******************************** */

template <typename T>
void TemporarilyProtectImmersedWalls3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;

    PLB_PRECONDITION(atomicBlocks.size() == 2);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int> *>(atomicBlocks[0]);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D *>(atomicBlocks[1]);
    PLB_ASSERT(flag);
    PLB_ASSERT(container);
    Dot3D location = flag->getLocation();
    ImmersedWallData3D<T> *wallData = dynamic_cast<ImmersedWallData3D<T> *>(container->getData());
    PLB_ASSERT(wallData);

    std::vector<Array<T, 3> > const &vertices = wallData->vertices;

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 3> const &vertex = vertices[i];
        Array<plint, 3> intPos(
            (plint)vertex[0] - location.x, (plint)vertex[1] - location.y,
            (plint)vertex[2] - location.z);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> zLim(
            (vertex[2] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        // x   x . x   x
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                for (plint dz = zLim[0]; dz <= zLim[1]; dz++) {
                    Array<plint, 3> pos(intPos + Array<plint, 3>(dx, dy, dz));
                    if (contained(pos[0], pos[1], pos[2], domain)) {
                        if (flag->get(pos[0], pos[1], pos[2]) == fluid) {
                            flag->get(pos[0], pos[1], pos[2]) = temporarilyProtect;
                        }
                    }
                }
            }
        }
    }
}

/* *************** Class RemoveProtectionFromImmersedWalls3D ******************************** */

template <typename T>
void RemoveProtectionFromImmersedWalls3D<T>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> atomicBlocks)
{
    using namespace freeSurfaceFlag;

    PLB_PRECONDITION(atomicBlocks.size() == 1);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int> *>(atomicBlocks[0]);
    PLB_ASSERT(flag);

    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                if (flag->get(iX, iY, iZ) == temporarilyProtect) {
                    flag->get(iX, iY, iZ) = fluid;
                }
            }
        }
    }
}

}  // namespace plb

#endif  // FREE_SURFACE_MODEL_3D_HH
