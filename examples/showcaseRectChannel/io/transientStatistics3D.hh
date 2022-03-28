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

#ifndef TRANSIENT_STATISTICS_3D_HH
#define TRANSIENT_STATISTICS_3D_HH

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "dataProcessors/dataAnalysisFunctional3D.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "io/imageWriter.h"
#include "io/plbFiles.h"
#include "io/transientStatistics3D.h"
#include "io/vtkDataOutput.h"
#include "multiBlock/multiBlock3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"

namespace plb {

/* ***************** Transient Statistics Manager ************************* */

template <typename T, template <typename U> class Descriptor>
TransientStatistics3D<T, Descriptor>::TransientStatistics3D(
    MultiBlockLattice3D<T, Descriptor> &lattice_, Box3D const &domain_) :
    lattice(lattice_)
{
#ifdef PLB_DEBUG
    bool intersectsWithSimulationDomain =
#endif
        intersect(domain_, lattice.getBoundingBox(), domain);
    PLB_ASSERT(intersectsWithSimulationDomain);
    enlargedDomain = domain;
    n = 0;
    isInitialized = 0;
    (void)memset(fieldIsRegistered, 0, sizeof fieldIsRegistered);
    (void)memset(fieldOperationIsRegistered, 0, sizeof fieldOperationIsRegistered);
    (void)memset(blocks, 0, sizeof blocks);
}

template <typename T, template <typename U> class Descriptor>
TransientStatistics3D<T, Descriptor>::TransientStatistics3D(
    TransientStatistics3D<T, Descriptor> const &rhs) :
    lattice(rhs.lattice),
    domain(rhs.domain),
    enlargedDomain(rhs.enlargedDomain),
    n(rhs.n),
    isInitialized(rhs.isInitialized)
{
    for (int iField = 0; iField < numFields; iField++) {
        fieldIsRegistered[iField] = rhs.fieldIsRegistered[iField];
    }

    for (int iField = 0; iField < numFields; iField++) {
        for (int iOperation = 0; iOperation < numOperations; iOperation++) {
            fieldOperationIsRegistered[iField][iOperation] =
                rhs.fieldOperationIsRegistered[iField][iOperation];
        }
    }

    for (int iField = 0; iField < numFields; iField++) {
        for (int iOperation = 0; iOperation < numOperations; iOperation++) {
            blocks[iField][iOperation] = 0;
            if (rhs.blocks[iField][iOperation] != 0) {
                blocks[iField][iOperation] = rhs.blocks[iField][iOperation]->clone();
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void TransientStatistics3D<T, Descriptor>::swap(TransientStatistics3D<T, Descriptor> &rhs)
{
    std::swap(lattice, rhs.lattice);
    std::swap(domain, rhs.domain);
    std::swap(enlargedDomain, rhs.enlargedDomain);
    std::swap(n, rhs.n);
    std::swap(isInitialized, rhs.isInitialized);

    for (int iField = 0; iField < numFields; iField++) {
        std::swap(fieldIsRegistered[iField], rhs.fieldIsRegistered[iField]);
    }

    for (int iField = 0; iField < numFields; iField++) {
        for (int iOperation = 0; iOperation < numOperations; iOperation++) {
            std::swap(
                fieldOperationIsRegistered[iField][iOperation],
                rhs.fieldOperationIsRegistered[iField][iOperation]);
        }
    }

    for (int iField = 0; iField < numFields; iField++) {
        for (int iOperation = 0; iOperation < numOperations; iOperation++) {
            std::swap(blocks[iField][iOperation], rhs.blocks[iField][iOperation]);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
TransientStatistics3D<T, Descriptor> &TransientStatistics3D<T, Descriptor>::operator=(
    TransientStatistics3D<T, Descriptor> const &rhs)
{
    TransientStatistics3D<T, Descriptor>(rhs).swap(*this);
    return *this;
}

template <typename T, template <typename U> class Descriptor>
TransientStatistics3D<T, Descriptor> *TransientStatistics3D<T, Descriptor>::clone() const
{
    return new TransientStatistics3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
TransientStatistics3D<T, Descriptor>::~TransientStatistics3D()
{
    for (int iField = 0; iField < numFields; iField++) {
        for (int iOperation = 0; iOperation < numOperations; iOperation++) {
            delete blocks[iField][iOperation];
        }
    }
}

template <typename T, template <typename U> class Descriptor>
bool TransientStatistics3D<T, Descriptor>::registerFieldOperation(
    std::string field, std::string operation)
{
    if (isInitialized) {  // No registering is allowed after initialization.
        return false;
    }

    int iField = fieldToId(field);
    PLB_ASSERT(iField >= 0);
    int iOperation = operationToId(operation);
    PLB_ASSERT(iOperation >= 0);

    fieldIsRegistered[iField] = 1;
    fieldOperationIsRegistered[iField][iOperation] = 1;
    if (iOperation == dev) {
        fieldOperationIsRegistered[iField][ave] =
            1;  // When the standard deviation is registered, averaging is registered also.
    }

    return true;
}

template <typename T, template <typename U> class Descriptor>
void TransientStatistics3D<T, Descriptor>::initialize()
{
    if (isInitialized) {
        return;
    }

    // We use an enlarged domain to compute the vorticity, whenever it is registered.
    if (fieldIsRegistered[vorticityX] || fieldIsRegistered[vorticityY]
        || fieldIsRegistered[vorticityZ] || fieldIsRegistered[vorticityNorm])
    {
#ifdef PLB_DEBUG
        bool intersectsWithSimulationDomain =
#endif
            intersect(domain.enlarge(1), lattice.getBoundingBox(), enlargedDomain);
        PLB_ASSERT(intersectsWithSimulationDomain);
    }

    for (int iField = 0; iField < numFields; iField++) {
        if (fieldIsRegistered[iField]) {  // There is at least one operation for this field.
            MultiScalarField3D<T> *field = computeField(iField);
            for (int iOperation = 0; iOperation < numOperations; iOperation++) {
                if (fieldOperationIsRegistered[iField][iOperation]) {
                    MultiScalarField3D<T> *fieldToRegister = 0;
                    if (iOperation == rms) {
                        fieldToRegister = computeAbsoluteValue(*field).release();
                    } else if (iOperation == dev) {
                        fieldToRegister = field->clone();
                        setToConstant(*fieldToRegister, fieldToRegister->getBoundingBox(), (T)0);
                    } else {
                        fieldToRegister = field->clone();
                    }

                    blocks[iField][iOperation] = fieldToRegister;
                }
            }
            delete field;
        }
    }

    n = 1;
    isInitialized = 1;
}

template <typename T, template <typename U> class Descriptor>
void TransientStatistics3D<T, Descriptor>::update()
{
    if (!isInitialized) {
        initialize();
        return;
    }

    n++;

    for (int iField = 0; iField < numFields; iField++) {
        if (fieldIsRegistered[iField]) {  // There is at least one operation for this field.
            MultiScalarField3D<T> *field = computeField(iField);
            for (int iOperation = 0; iOperation < numOperations; iOperation++) {
                if (fieldOperationIsRegistered[iField][iOperation]) {
                    std::vector<MultiBlock3D *> args;
                    args.push_back(field);
                    if (iOperation == dev) {
                        args.push_back(blocks[iField][ave]);  // ave must be less than dev, so that
                                                              // the average is already updated.
                    }
                    args.push_back(blocks[iField][iOperation]);

                    BoxProcessingFunctional3D *functional = 0;
                    switch (iOperation) {
                    case min:
                        functional = new UpdateMinScalarTransientStatistics3D<T>();
                        break;
                    case max:
                        functional = new UpdateMaxScalarTransientStatistics3D<T>();
                        break;
                    case ave:
                        functional = new UpdateAveScalarTransientStatistics3D<T>(n);
                        break;
                    case rms:
                        functional = new UpdateRmsScalarTransientStatistics3D<T>(n);
                        break;
                    case dev:
                        functional = new UpdateDevScalarTransientStatistics3D<T>(n);
                        break;
                    default:
                        functional = 0;
                        break;
                    }

                    applyProcessingFunctional(functional, field->getBoundingBox(), args);
                }
            }
            delete field;
        }
    }
}

template <typename T, template <typename U> class Descriptor>
MultiScalarField3D<T> &TransientStatistics3D<T, Descriptor>::get(
    std::string field, std::string operation) const
{
    int iField = fieldToId(field);
    PLB_ASSERT(iField >= 0);
    int iOperation = operationToId(operation);
    PLB_ASSERT(iOperation >= 0);

    return *blocks[iField][iOperation];
}

template <typename T, template <typename U> class Descriptor>
void TransientStatistics3D<T, Descriptor>::output(
    std::string path, std::string domainName, plint iteration, plint namePadding, T dx, T dt,
    Array<T, 3> const &physicalLocation, T rho, T pressureOffset, T rhoLB)
{
    if (!isInitialized) {
        initialize();
    }

    for (int iField = 0; iField < numFields; iField++) {
        for (int iOperation = 0; iOperation < numOperations; iOperation++) {
            if (fieldOperationIsRegistered[iField][iOperation]) {
                std::string fileName =
                    getFileName(path, iField, iOperation, domainName, iteration, namePadding);
                VtkImageOutput3D<T> vtkOut(fileName, dx, physicalLocation);
                std::string field = idToField(iField);
                T scalingFactor = getScalingFactor(iField, dx, dt, rho);
                T offset = getOffset(iField, dx, dt, rho, pressureOffset, rhoLB);
                vtkOut.writeData(*blocks[iField][iOperation], field, scalingFactor, offset);
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void TransientStatistics3D<T, Descriptor>::saveState(
    plint iteration, FileName xmlFileName, FileName baseFileName, plint fileNamePadding)
{
    if (!isInitialized) {
        initialize();
    }

    std::string fname_base = createFileName(baseFileName.get() + "_", iteration, fileNamePadding);

    bool saveDynamicContent = false;
    for (int iField = 0; iField < numFields; iField++) {
        if (fieldIsRegistered[iField]) {
            std::string field = idToField(iField);
            for (int iOperation = 0; iOperation < numOperations; iOperation++) {
                if (fieldOperationIsRegistered[iField][iOperation]) {
                    std::string operation = idToOperation(iOperation);
                    std::string fname = fname_base + "_" + field + "_" + operation;
                    parallelIO::save(*blocks[iField][iOperation], fname, saveDynamicContent);
                }
            }
        }
    }

    XMLwriter restart;
    XMLwriter &entry = restart["continue"]["transientStatistics"];
    entry["name"].setString(FileName(fname_base).defaultPath(global::directories().getOutputDir()));
    entry["iteration"].set(iteration);
    entry["n"].set(n);
    restart.print(xmlFileName);
}

template <typename T, template <typename U> class Descriptor>
void TransientStatistics3D<T, Descriptor>::loadState(plint &iteration, FileName xmlFileName)
{
    XMLreader restart(xmlFileName.get());
    std::string fname_base;
    restart["continue"]["transientStatistics"]["name"].read(fname_base);
    restart["continue"]["transientStatistics"]["iteration"].read(iteration);
    restart["continue"]["transientStatistics"]["n"].read(n);

    bool saveDynamicContent = false;
    for (int iField = 0; iField < numFields; iField++) {
        if (fieldIsRegistered[iField]) {
            std::string field = idToField(iField);
            for (int iOperation = 0; iOperation < numOperations; iOperation++) {
                if (fieldOperationIsRegistered[iField][iOperation]) {
                    std::string operation = idToOperation(iOperation);
                    std::string fname = fname_base + "_" + field + "_" + operation;
                    parallelIO::load(fname, *blocks[iField][iOperation], saveDynamicContent);
                }
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
int TransientStatistics3D<T, Descriptor>::fieldToId(std::string field) const
{
    if (field == "velocityX") {
        return velocityX;
    } else if (field == "velocityY") {
        return velocityY;
    } else if (field == "velocityZ") {
        return velocityZ;
    } else if (field == "velocityNorm") {
        return velocityNorm;
    } else if (field == "density") {
        return density;
    } else if (field == "pressure") {
        return pressure;
    } else if (field == "vorticityX") {
        return vorticityX;
    } else if (field == "vorticityY") {
        return vorticityY;
    } else if (field == "vorticityZ") {
        return vorticityZ;
    } else if (field == "vorticityNorm") {
        return vorticityNorm;
    } else {
        return -1;
    }
}

template <typename T, template <typename U> class Descriptor>
int TransientStatistics3D<T, Descriptor>::operationToId(std::string operation) const
{
    if (operation == "min") {
        return min;
    } else if (operation == "max") {
        return max;
    } else if (operation == "ave") {
        return ave;
    } else if (operation == "rms") {
        return rms;
    } else if (operation == "dev") {
        return dev;
    } else {
        return -1;
    }
}

template <typename T, template <typename U> class Descriptor>
std::string TransientStatistics3D<T, Descriptor>::idToField(int iField) const
{
    PLB_ASSERT(iField >= 0);

    switch (iField) {
    case velocityX:
        return std::string("velocityX");
    case velocityY:
        return std::string("velocityY");
    case velocityZ:
        return std::string("velocityZ");
    case velocityNorm:
        return std::string("velocityNorm");
    case density:
        return std::string("density");
    case pressure:
        return std::string("pressure");
    case vorticityX:
        return std::string("vorticityX");
    case vorticityY:
        return std::string("vorticityY");
    case vorticityZ:
        return std::string("vorticityZ");
    case vorticityNorm:
        return std::string("vorticityNorm");
    default:
        return std::string("error");
    }
}

template <typename T, template <typename U> class Descriptor>
std::string TransientStatistics3D<T, Descriptor>::idToOperation(int iOperation) const
{
    PLB_ASSERT(iOperation >= 0);

    switch (iOperation) {
    case min:
        return std::string("min");
    case max:
        return std::string("max");
    case ave:
        return std::string("ave");
    case rms:
        return std::string("rms");
    case dev:
        return std::string("dev");
    default:
        return std::string("error");
    }
}

template <typename T, template <typename U> class Descriptor>
T TransientStatistics3D<T, Descriptor>::getScalingFactor(int iField, T dx, T dt, T rho) const
{
    PLB_ASSERT(iField >= 0);

    T scalingFactor = 0.0;

    switch (iField) {
    case velocityX:
    case velocityY:
    case velocityZ:
    case velocityNorm:
        scalingFactor = dx / dt;
        break;
    case density:
        scalingFactor = rho;
        break;
    case pressure:
        scalingFactor = rho * dx * dx / (dt * dt) * Descriptor<T>::cs2;
        break;
    case vorticityX:
    case vorticityY:
    case vorticityZ:
    case vorticityNorm:
        scalingFactor = 1.0 / dt;
        break;
    default:
        scalingFactor = 0.0;
        break;
    }

    return scalingFactor;
}

template <typename T, template <typename U> class Descriptor>
T TransientStatistics3D<T, Descriptor>::getOffset(
    int iField, T dx, T dt, T rho, T pressureOffset, T rhoLB) const
{
    PLB_ASSERT(iField >= 0);

    T offset = 0.0;

    switch (iField) {
    case pressure:
        offset = pressureOffset - rho * dx * dx / (dt * dt) * Descriptor<T>::cs2 * rhoLB;
        break;
    default:
        offset = 0.0;
        break;
    }

    return offset;
}

template <typename T, template <typename U> class Descriptor>
MultiScalarField3D<T> *TransientStatistics3D<T, Descriptor>::computeField(int iField) const
{
    PLB_ASSERT(iField >= 0);

    MultiScalarField3D<T> *field = 0;

    switch (iField) {
    case velocityX:
        field = computeVelocityComponent(lattice, domain, 0).release();
        break;
    case velocityY:
        field = computeVelocityComponent(lattice, domain, 1).release();
        break;
    case velocityZ:
        field = computeVelocityComponent(lattice, domain, 2).release();
        break;
    case velocityNorm:
        field = computeVelocityNorm(lattice, domain).release();
        break;
    case density:
    case pressure:
        field = computeDensity(lattice, domain).release();
        break;
    case vorticityX:
        field = extractComponent(
                    *computeVorticity(*computeVelocity(lattice, enlargedDomain)), domain, 0)
                    .release();
        break;
    case vorticityY:
        field = extractComponent(
                    *computeVorticity(*computeVelocity(lattice, enlargedDomain)), domain, 1)
                    .release();
        break;
    case vorticityZ:
        field = extractComponent(
                    *computeVorticity(*computeVelocity(lattice, enlargedDomain)), domain, 2)
                    .release();
        break;
    case vorticityNorm:
        field = computeNorm(*computeVorticity(*computeVelocity(lattice, enlargedDomain)), domain)
                    .release();
        break;
    default:
        field = 0;
        break;
    }

    return field;
}

template <typename T, template <typename U> class Descriptor>
std::string TransientStatistics3D<T, Descriptor>::getFileName(
    std::string path, int iField, int iOperation, std::string domainName, plint iteration,
    plint namePadding) const
{
    std::string field = idToField(iField);
    std::string operation = idToOperation(iOperation);
    std::string name =
        createFileName(field + "_" + operation + "_" + domainName + "_", iteration, namePadding);
    FileName fileName;
    fileName.setPath(path);
    fileName.setName(name);
    return fileName.get();
}

}  // namespace plb

#endif  // TRANSIENT_STATISTICS_3D_HH
