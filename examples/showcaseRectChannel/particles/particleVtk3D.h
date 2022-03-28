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

#ifndef PARTICLE_VTK_3D_H
#define PARTICLE_VTK_3D_H

#include <string>
#include <vector>

#include "core/functions.h"
#include "core/globalDefs.h"
#include "offLattice/triangleBoundary3D.h"
#include "particles/multiParticleField3D.h"

namespace plb {

template <typename T, template <typename U> class Descriptor, class BoundaryType>
void writeSurfaceVTK(
    TriangleBoundary3D<T> const &boundary,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles,
    std::vector<std::string> const &scalars, std::vector<std::string> const &vectors,
    std::string const &fName, bool dynamicMesh, plint tag);

template <typename T, template <typename U> class Descriptor, class BoundaryType>
void writeSurfaceVTK(
    TriangleBoundary3D<T> const &boundary,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles,
    std::vector<std::string> const &scalars, std::vector<std::string> const &vectors,
    std::string const &fName, bool dynamicMesh, plint tag, std::vector<T> const &scalarFactor,
    std::vector<T> const &vectorFactor);

template <typename T, template <typename U> class Descriptor, class BoundaryType>
void vtkForVertices(
    std::vector<Particle3D<T, Descriptor> *> const &particles,
    TriangleBoundary3D<T> const &boundary, std::vector<std::string> const &scalars,
    std::vector<std::string> const &vectors, std::string fName, bool dynamicMesh, plint tag,
    std::vector<T> const &scalarFactor, std::vector<T> const &vectorFactor);

template <typename T, template <typename U> class Descriptor, class BoundaryType>
void writeVertexAsciiData(
    TriangleBoundary3D<T> const &boundary,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles,
    std::vector<std::string> const &scalars, std::vector<std::string> const &vectors,
    std::string const &fName, bool dynamicMesh, bool printHeader,
    std::vector<T> const &scalarFactor, std::vector<T> const &vectorFactor);

template <typename T, template <typename U> class Descriptor, class BoundaryType>
void vertexAsciiData(
    std::vector<Particle3D<T, Descriptor> *> const &particles,
    TriangleBoundary3D<T> const &boundary, std::vector<std::string> const &scalars,
    std::vector<std::string> const &vectors, std::string fName, bool dynamicMesh, bool printHeader,
    std::vector<T> const &scalarFactor, std::vector<T> const &vectorFactor);

template <typename T, template <typename U> class Descriptor>
void writeAsciiParticlePos(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    T deltaX, Array<T, 3> const &offset);

template <typename T, template <typename U> class Descriptor>
void writeAsciiParticlePos(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    T deltaX = T(1))
{
    Array<T, 3> offset;
    offset.resetToZero();
    writeAsciiParticlePos(particles, fName, deltaX, offset);
}

template <typename T, template <typename U> class Descriptor>
void writeParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    T deltaX, Array<T, 3> const &offset, pluint maxNumParticlesToWrite = 0);

template <typename T, template <typename U> class Descriptor>
void writeParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    T deltaX = T(1), pluint maxNumParticlesToWrite = 0)
{
    Array<T, 3> offset;
    offset.resetToZero();
    writeParticleVtk(particles, fName, deltaX, offset, maxNumParticlesToWrite);
}

template <typename T, template <typename U> class Descriptor>
void writeParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    std::map<plint, std::string> const &additionalScalars,
    std::map<plint, std::string> const &additionalVectors, T deltaX, Array<T, 3> const &offset,
    pluint maxNumParticlesToWrite = 0);

template <typename T, template <typename U> class Descriptor>
void writeParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    std::map<plint, std::string> const &additionalScalars,
    std::map<plint, std::string> const &additionalVectors, T deltaX = T(1),
    pluint maxNumParticlesToWrite = 0)
{
    Array<T, 3> offset;
    offset.resetToZero();
    writeParticleVtk(
        particles, fName, additionalScalars, additionalVectors, deltaX, offset,
        maxNumParticlesToWrite);
}

template <typename T, template <typename U> class Descriptor>
void writeSelectedParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    Box3D const &domain, util::SelectInt const &tags, T deltaX, Array<T, 3> const &offset);

template <typename T, template <typename U> class Descriptor>
void writeSelectedParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    Box3D const &domain, util::SelectInt const &tags, T deltaX = T(1))
{
    Array<T, 3> offset;
    offset.resetToZero();
    writeSelectedParticleVtk(particles, fName, domain, tags, deltaX, offset);
}

template <typename T, template <typename U> class Descriptor>
void writeSelectedParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    Box3D const &domain, plint tag, T deltaX, Array<T, 3> const &offset)
{
    writeSelectedParticleVtk(particles, fName, domain, util::SelectConstInt(tag), deltaX, offset);
}

template <typename T, template <typename U> class Descriptor>
void writeSelectedParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    Box3D const &domain, plint tag, T deltaX = T(1))
{
    Array<T, 3> offset;
    offset.resetToZero();
    writeSelectedParticleVtk(particles, fName, domain, tag, deltaX, offset);
}

template <typename T, template <typename U> class Descriptor>
void writeSelectedParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    Box3D const &domain, util::SelectInt const &tags,
    std::map<plint, std::string> const &additionalScalars,
    std::map<plint, std::string> const &additionalVectors, T deltaX, Array<T, 3> const &offset);

template <typename T, template <typename U> class Descriptor>
void writeSelectedParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    Box3D const &domain, util::SelectInt const &tags,
    std::map<plint, std::string> const &additionalScalars,
    std::map<plint, std::string> const &additionalVectors, T deltaX = T(1))
{
    Array<T, 3> offset;
    offset.resetToZero();
    writeSelectedParticleVtk(
        particles, fName, domain, tags, additionalScalars, additionalVectors, deltaX, offset);
}

template <typename T, template <typename U> class Descriptor>
void writeSelectedParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    Box3D const &domain, plint tag, std::map<plint, std::string> const &additionalScalars,
    std::map<plint, std::string> const &additionalVectors, T deltaX, Array<T, 3> const &offset)
{
    writeSelectedParticleVtk(
        particles, fName, domain, util::SelectConstInt(tag), additionalScalars, additionalVectors,
        deltaX, offset);
}

template <typename T, template <typename U> class Descriptor>
void writeSelectedParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    Box3D const &domain, plint tag, std::map<plint, std::string> const &additionalScalars,
    std::map<plint, std::string> const &additionalVectors, T deltaX = T(1))
{
    Array<T, 3> offset;
    offset.resetToZero();
    writeSelectedParticleVtk(
        particles, fName, domain, tag, additionalScalars, additionalVectors, deltaX, offset);
}

}  // namespace plb

#endif  // PARTICLE_VTK_3D_H
