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

#ifndef TRANSIENT_STATISTICS_3D_H
#define TRANSIENT_STATISTICS_3D_H

#include <string>

#include "core/array.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "io/plbFiles.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

/* ***************** Transient Statistics Manager ************************* */

template <typename T, template <typename U> class Descriptor>
class TransientStatistics3D {
public:
    TransientStatistics3D(MultiBlockLattice3D<T, Descriptor> &lattice_, Box3D const &domain_);
    TransientStatistics3D(TransientStatistics3D<T, Descriptor> const &rhs);
    void swap(TransientStatistics3D<T, Descriptor> &rhs);
    TransientStatistics3D<T, Descriptor> &operator=(
        TransientStatistics3D<T, Descriptor> const &rhs);
    TransientStatistics3D<T, Descriptor> *clone() const;
    ~TransientStatistics3D();
    // Field must be one of:
    //
    // "velocityX", "velocityY", "velocityZ", "velocityNorm", "density", "pressure",
    //     "vorticityX", "vorticityY", "vorticityZ", "vorticityNorm"
    //
    // The "density" entry above is mainly to treat advection-diffusion equations.
    //
    // Operation must be one of:
    // "min", "max", "ave" (for mean value), "rms", "dev" (for standard deviation).
    bool registerFieldOperation(std::string field, std::string operation);
    void initialize();
    void update();
    MultiScalarField3D<T> &get(std::string field, std::string operation) const;
    // "rho" is the actual fluid density in physical units (or the scaling factor for an
    // advection-diffusion equation). "pressureOffset" is the ambient pressure in physical units.
    // "rhoLB" is the lattice-Boltzmann density (1 by default).
    void output(
        std::string path, std::string domainName, plint iteration, plint namePadding, T dx, T dt,
        Array<T, 3> const &physicalLocation, T rho = T(1), T pressureOffset = T(0), T rhoLB = T(1));
    // Save the current state of all the registered transient statistics for restarting.
    void saveState(
        plint iteration, FileName xmlFileName, FileName baseFileName, plint fileNamePadding = 8);
    // Save the current state of all the registered transient statistics for restarting.
    void loadState(plint &iteration, FileName xmlFileName);

private:
    int fieldToId(std::string field) const;
    int operationToId(std::string operation) const;
    std::string idToField(int iField) const;
    std::string idToOperation(int iOperation) const;
    T getScalingFactor(int iField, T dx, T dt, T rho) const;
    T getOffset(int iField, T dx, T dt, T rho, T pressureOffset, T rhoLB) const;
    MultiScalarField3D<T> *computeField(int iField) const;
    std::string getFileName(
        std::string path, int iField, int iOperation, std::string domainName, plint iteration,
        plint namePadding) const;

private:
    enum { numFields = 10 };
    enum {
        velocityX,
        velocityY,
        velocityZ,
        velocityNorm,
        density,
        pressure,
        vorticityX,
        vorticityY,
        vorticityZ,
        vorticityNorm
    };
    enum { numOperations = 5 };
    enum { min, max, ave, rms, dev };

private:
    MultiBlockLattice3D<T, Descriptor> &lattice;  // Reference to the lattice of the simulation.
    Box3D domain;                      // Domain to compute and output transient statistics.
    Box3D enlargedDomain;              // Domain to compute vorticity.
    plint n;                           // Number of times the statistics are updated.
    int isInitialized;                 // Is the data structure initialized or not.
    int fieldIsRegistered[numFields];  // Array of all registered fields.
    int fieldOperationIsRegistered[numFields][numOperations];  // Table of all registered fields and
                                                               // operations.
    MultiScalarField3D<T> *blocks[numFields][numOperations];   // All scalar fields to operate on.
};

}  // namespace plb

#endif  // TRANSIENT_STATISTICS_3D_H
