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

#ifndef PARTICLE_VTK_3D_HH
#define PARTICLE_VTK_3D_HH

#include <cstdio>
#include <cstdlib>
#include <set>

#include "core/globalDefs.h"
#include "particles/particleNonLocalTransfer3D.h"
#include "particles/particleVtk3D.h"
#include "sitmo/prng_engine.hpp"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void writeSurfaceVTK(
    TriangleBoundary3D<T> const &boundary,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles,
    std::vector<std::string> const &scalars, std::vector<std::string> const &vectors,
    std::string const &fName, bool dynamicMesh, plint tag)
{
    writeSurfaceVTK(
        boundary, particles, scalars, vectors, fName, dynamicMesh, tag, std::vector<T>(),
        std::vector<T>());
}

template <typename T, template <typename U> class Descriptor>
void writeSurfaceVTK(
    TriangleBoundary3D<T> const &boundary,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles,
    std::vector<std::string> const &scalars, std::vector<std::string> const &vectors,
    std::string const &fName, bool dynamicMesh, plint tag, std::vector<T> const &scalarFactor,
    std::vector<T> const &vectorFactor)
{
    writeSurfaceVTK(
        boundary, particles, particles.getBoundingBox(), scalars, vectors, fName, dynamicMesh, tag,
        scalarFactor, vectorFactor);
}

template <typename T, template <typename U> class Descriptor>
void writeSurfaceVTK(
    TriangleBoundary3D<T> const &boundary,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, Box3D const &bb,
    std::vector<std::string> const &scalars, std::vector<std::string> const &vectors,
    std::string const &fName, bool dynamicMesh, plint tag, std::vector<T> const &scalarFactor,
    std::vector<T> const &vectorFactor)
{
    SparseBlockStructure3D blockStructure(bb);
    blockStructure.addBlock(bb, 0);
    plint envelopeWidth = 1;
    MultiBlockManagement3D serialMultiBlockManagement(
        blockStructure, new OneToOneThreadAttribution, envelopeWidth);

    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > multiSerialParticles(
        serialMultiBlockManagement, defaultMultiBlockPolicy3D().getCombinedStatistics());

    copy(particles, bb, multiSerialParticles, bb);
    if (global::mpi().isMainProcessor()) {
        ParticleField3D<T, Descriptor> &atomicSerialParticles =
            dynamic_cast<ParticleField3D<T, Descriptor> &>(multiSerialParticles.getComponent(0));

        std::vector<Particle3D<T, Descriptor> *> found;
        SmartBulk3D oneBlockBulk(serialMultiBlockManagement, 0);
        atomicSerialParticles.findParticles(oneBlockBulk.toLocal(bb), found);

        vtkForVertices(
            found, boundary, scalars, vectors, fName, dynamicMesh, tag, scalarFactor, vectorFactor);
    }
}

template <typename T, template <typename U> class Descriptor>
void vtkForVertices(
    std::vector<Particle3D<T, Descriptor> *> const &particles,
    TriangleBoundary3D<T> const &boundary, std::vector<std::string> const &scalars,
    std::vector<std::string> const &vectors, std::string fName, bool dynamicMesh, plint tag,
    std::vector<T> const &scalarFactor, std::vector<T> const &vectorFactor)
{
    PLB_ASSERT(scalarFactor.empty() || scalarFactor.size() == scalars.size());
    PLB_ASSERT(vectorFactor.empty() || vectorFactor.size() == vectors.size());
    if (dynamicMesh) {
        boundary.pushSelect(0, 1);  // 0=Open, 1=Dynamic.
    } else {
        boundary.pushSelect(0, 0);  // Open, Static.
    }
    TriangularSurfaceMesh<T> const &mesh = boundary.getMesh();
    // If this assertion fails, a likely explanation is that the margin of your sparse block
    // structure is too small, and one of the particles was outside the allocated domain.
    // PLB_PRECONDITION((plint)particles.size() == mesh.getNumVertices());
    if ((plint)particles.size() != mesh.getNumVertices()) {
        pcout << "Warning: in vtkForVertices, the numer of particles doesn't match the number"
              << std::endl;
        pcout << "of vertices. There might be black spots in the produced VTK file." << std::endl;
    }

    std::ofstream ofile(fName.c_str());
    ofile.precision(10);
    std::scientific(ofile);
    ofile << "# vtk DataFile Version 3.0\n";
    ofile << "Surface mesh created with Palabos\n";
    ofile << "ASCII\n";
    ofile << "DATASET UNSTRUCTURED_GRID\n";

    ofile << "POINTS " << particles.size() << (sizeof(T) == sizeof(double) ? " double" : " float")
          << "\n";

    std::vector<std::vector<T> > scalarData(scalars.size());
    for (pluint iScalar = 0; iScalar < scalars.size(); ++iScalar) {
        scalarData[iScalar].resize(mesh.getNumVertices());
    }
    std::vector<std::vector<Array<T, 3> > > vectorData(vectors.size());
    for (pluint iVector = 0; iVector < vectors.size(); ++iVector) {
        vectorData[iVector].resize(mesh.getNumVertices());
    }
    std::vector<Array<T, 3> > posVect(mesh.getNumVertices());
    // Default initialize all data, just in case every vertex hasn't
    // an associated particle.
    for (plint i = 0; i < mesh.getNumVertices(); ++i) {
        for (pluint iScalar = 0; iScalar < scalars.size(); ++iScalar) {
            scalarData[iScalar][i] = T();
        }
        for (pluint iVector = 0; iVector < vectors.size(); ++iVector) {
            vectorData[iVector][i].resetToZero();
        }
        posVect[i].resetToZero();
    }
    for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
        plint iVertex = particles[iParticle]->getTag();
        posVect[iVertex] = particles[iParticle]->getPosition();
        posVect[iVertex] *= boundary.getDx();
        posVect[iVertex] += boundary.getPhysicalLocation();
        for (pluint iScalar = 0; iScalar < scalars.size(); ++iScalar) {
            T scalar;
            particles[iParticle]->getScalar(iScalar, scalar);
            if (!scalarFactor.empty()) {
                scalar *= scalarFactor[iScalar];
            }
            scalarData[iScalar][iVertex] = scalar;
        }
        for (pluint iVector = 0; iVector < vectors.size(); ++iVector) {
            Array<T, 3> vector;
            particles[iParticle]->getVector(iVector, vector);
            if (!vectorFactor.empty()) {
                vector *= vectorFactor[iVector];
            }
            vectorData[iVector][iVertex] = vector;
        }
    }

    for (plint iVertex = 0; iVertex < mesh.getNumVertices(); ++iVertex) {
        ofile << posVect[iVertex][0] << " " << posVect[iVertex][1] << " " << posVect[iVertex][2]
              << "\n";
    }
    ofile << "\n";

    plint numWallTriangles = 0;
    for (plint iTriangle = 0; iTriangle < mesh.getNumTriangles(); ++iTriangle) {
        if (tag < 0 || boundary.getTag(iTriangle) == tag) {
            ++numWallTriangles;
        }
    }

    ofile << "CELLS " << numWallTriangles << " " << 4 * numWallTriangles << "\n";

    for (plint iTriangle = 0; iTriangle < mesh.getNumTriangles(); ++iTriangle) {
        plint i0 = mesh.getVertexId(iTriangle, 0);
        plint i1 = mesh.getVertexId(iTriangle, 1);
        plint i2 = mesh.getVertexId(iTriangle, 2);
        if (tag < 0 || boundary.getTag(iTriangle) == tag) {
            ofile << "3 " << i0 << " " << i1 << " " << i2 << "\n";
        }
    }
    ofile << "\n";

    ofile << "CELL_TYPES " << numWallTriangles << "\n";
    for (plint i = 0; i < numWallTriangles; ++i) {
        ofile << "5\n";
    }
    ofile << "\n";

    ofile << "POINT_DATA " << mesh.getNumVertices() << "\n";
    for (pluint iVector = 0; iVector < vectors.size(); ++iVector) {
        ofile << "VECTORS " << vectors[iVector]
              << (sizeof(T) == sizeof(double) ? " double" : " float") << "\n";
        for (plint iVertex = 0; iVertex < mesh.getNumVertices(); ++iVertex) {
            ofile << vectorData[iVector][iVertex][0] << " " << vectorData[iVector][iVertex][1]
                  << " " << vectorData[iVector][iVertex][2] << "\n";
        }
        ofile << "\n";
    }

    for (pluint iScalar = 0; iScalar < scalars.size(); ++iScalar) {
        ofile << "SCALARS " << scalars[iScalar]
              << (sizeof(T) == sizeof(double) ? " double" : " float") << " 1\n"
              << "LOOKUP_TABLE default\n";
        for (plint iVertex = 0; iVertex < mesh.getNumVertices(); ++iVertex) {
            ofile << scalarData[iScalar][iVertex] << "\n";
        }
        ofile << "\n";
    }
    // Restore mesh selection which was active before calling
    //   this function.
    boundary.popSelect();
}

template <typename T, template <typename U> class Descriptor>
void writeVertexAsciiData(
    TriangleBoundary3D<T> const &boundary,
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles,
    std::vector<std::string> const &scalars, std::vector<std::string> const &vectors,
    std::string const &fName, bool dynamicMesh, bool printHeader,
    std::vector<T> const &scalarFactor, std::vector<T> const &vectorFactor)
{
    SparseBlockStructure3D blockStructure(particles.getBoundingBox());
    blockStructure.addBlock(particles.getBoundingBox(), 0);
    plint envelopeWidth = 1;
    MultiBlockManagement3D serialMultiBlockManagement(
        blockStructure, new OneToOneThreadAttribution, envelopeWidth);

    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > multiSerialParticles(
        serialMultiBlockManagement, defaultMultiBlockPolicy3D().getCombinedStatistics());

    copy(particles, particles.getBoundingBox(), multiSerialParticles, particles.getBoundingBox());
    if (global::mpi().isMainProcessor()) {
        ParticleField3D<T, Descriptor> &atomicSerialParticles =
            dynamic_cast<ParticleField3D<T, Descriptor> &>(multiSerialParticles.getComponent(0));

        std::vector<Particle3D<T, Descriptor> *> found;
        SmartBulk3D oneBlockBulk(serialMultiBlockManagement, 0);
        atomicSerialParticles.findParticles(
            oneBlockBulk.toLocal(particles.getBoundingBox()), found);

        vertexAsciiData(
            found, boundary, scalars, vectors, fName, dynamicMesh, printHeader, scalarFactor,
            vectorFactor);
    }
}

template <typename T, template <typename U> class Descriptor>
void vertexAsciiData(
    std::vector<Particle3D<T, Descriptor> *> const &particles,
    TriangleBoundary3D<T> const &boundary, std::vector<std::string> const &scalars,
    std::vector<std::string> const &vectors, std::string fName, bool dynamicMesh, bool printHeader,
    std::vector<T> const &scalarFactor, std::vector<T> const &vectorFactor)
{
    PLB_ASSERT(scalarFactor.empty() || scalarFactor.size() == scalars.size());
    PLB_ASSERT(vectorFactor.empty() || vectorFactor.size() == vectors.size());
    if (dynamicMesh) {
        boundary.pushSelect(0, 1);  // 0=Open, 1=Dynamic.
    } else {
        boundary.pushSelect(0, 0);  // Open, Static.
    }
    TriangularSurfaceMesh<T> const &mesh = boundary.getMesh();
    PLB_PRECONDITION((plint)particles.size() == mesh.getNumVertices());
    std::ofstream ofile(fName.c_str());
    ofile.precision(10);
    std::scientific(ofile);
    if (printHeader) {
        ofile << "Pos[0] Pos[1] Pos[2]";
        for (pluint iScalar = 0; iScalar < scalars.size(); ++iScalar) {
            ofile << " " << scalars[iScalar];
        }
        for (pluint iVector = 0; iVector < vectors.size(); ++iVector) {
            ofile << " " << vectors[iVector] << "[0]";
            ofile << " " << vectors[iVector] << "[1]";
            ofile << " " << vectors[iVector] << "[2]";
        }
        ofile << "\n";
    }
    for (pluint iParticle = 0; iParticle < particles.size(); ++iParticle) {
        Array<T, 3> position = particles[iParticle]->getPosition();
        ofile << position[0] << " " << position[1] << " " << position[2];
        for (pluint iScalar = 0; iScalar < scalars.size(); ++iScalar) {
            T scalar;
            particles[iParticle]->getScalar(iScalar, scalar);
            if (!scalarFactor.empty()) {
                scalar *= scalarFactor[iScalar];
            }
            ofile << " " << scalar;
        }
        for (pluint iVector = 0; iVector < vectors.size(); ++iVector) {
            Array<T, 3> vector;
            particles[iParticle]->getVector(iVector, vector);
            if (!vectorFactor.empty()) {
                vector *= vectorFactor[iVector];
            }
            ofile << " " << vector[0] << " " << vector[1] << " " << vector[2];
        }
        ofile << "\n";
    }
    boundary.popSelect();
}

template <typename T, template <typename U> class Descriptor>
void writeAsciiParticlePos(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    T deltaX, Array<T, 3> const &offset)
{
    SparseBlockStructure3D blockStructure(particles.getBoundingBox());
    blockStructure.addBlock(particles.getBoundingBox(), 0);
    plint envelopeWidth = 1;
    MultiBlockManagement3D serialMultiBlockManagement(
        blockStructure, new OneToOneThreadAttribution, envelopeWidth);

    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > multiSerialParticles(
        serialMultiBlockManagement, defaultMultiBlockPolicy3D().getCombinedStatistics());

    copy(particles, particles.getBoundingBox(), multiSerialParticles, particles.getBoundingBox());
    if (global::mpi().isMainProcessor()) {
        std::ofstream ofile(fName.c_str());
        ofile.precision(10);
        std::scientific(ofile);
        ParticleField3D<T, Descriptor> &atomicSerialParticles =
            dynamic_cast<ParticleField3D<T, Descriptor> &>(multiSerialParticles.getComponent(0));

        std::vector<Particle3D<T, Descriptor> *> found;
        SmartBulk3D oneBlockBulk(serialMultiBlockManagement, 0);
        atomicSerialParticles.findParticles(
            oneBlockBulk.toLocal(particles.getBoundingBox()), found);
        for (pluint iParticle = 0; iParticle < found.size(); ++iParticle) {
            Array<T, 3> pos(found[iParticle]->getPosition());
            pos = deltaX * pos + offset;
            ofile << pos[0] << "," << pos[1] << "," << pos[2] << "\n";
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void particleVtkSerialImplementation(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &multiSerialParticles,
    std::string const &fName, std::map<plint, std::string> const &additionalScalars,
    std::map<plint, std::string> const &additionalVectors, T deltaX, Array<T, 3> const &offset,
    pluint maxNumParticlesToWrite)
{
    if (global::mpi().isMainProcessor()) {
        ParticleField3D<T, Descriptor> &atomicSerialParticles =
            dynamic_cast<ParticleField3D<T, Descriptor> &>(multiSerialParticles.getComponent(0));

        std::vector<Particle3D<T, Descriptor> *> found;
        SmartBulk3D oneBlockBulk(multiSerialParticles.getMultiBlockManagement(), 0);
        atomicSerialParticles.findParticles(
            oneBlockBulk.toLocal(multiSerialParticles.getBoundingBox()), found);
        pluint numParticles = found.size();
        if (numParticles == 0) {
            return;
        }
        pluint numParticlesToWrite =
            maxNumParticlesToWrite == 0
                ? numParticles
                : (maxNumParticlesToWrite > numParticles ? numParticles : maxNumParticlesToWrite);

        std::set<pluint> iParticles;
        if (numParticlesToWrite != numParticles) {
            sitmo::prng_engine eng;
            pluint count = 0;
            pluint maxIter = 1000000000;
            pluint iter = 0;
            while (count < numParticlesToWrite && iter < maxIter) {
                T randomValue = (T)eng() / ((T)sitmo::prng_engine::max() + 1.0);
                auto ret = iParticles.insert((pluint)(randomValue * (numParticles - 1)));
                if (ret.second) {
                    count++;
                }
                iter++;
            }
            if (count < numParticlesToWrite) {
                std::cout << "particleVtkSerialImplementation(): failed to randomly select "
                          << numParticlesToWrite << " particles after " << maxIter
                          << " iterations. Writing " << count << " particles instead." << std::endl;
            }
        }

        FILE *fp = fopen(fName.c_str(), "w");
        PLB_ASSERT(fp != 0);
        fprintf(fp, "# vtk DataFile Version 3.0\n");
        fprintf(fp, "Particle file created by Palabos\n");
        fprintf(fp, "ASCII\n");
        fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
        unsigned long ulNumParticlesToWrite = (unsigned long)numParticlesToWrite;
        fprintf(fp, "POINTS %12lu double\n", ulNumParticlesToWrite);
        if (numParticlesToWrite == numParticles) {
            for (pluint iParticle = 0; iParticle < numParticles; iParticle++) {
                Array<T, 3> pos(found[iParticle]->getPosition());
                pos = deltaX * pos + offset;
                Array<double, 3> doublePos(pos);
                fprintf(fp, "% .10e % .10e % .10e\n", doublePos[0], doublePos[1], doublePos[2]);
            }
        } else {
            for (pluint iParticle : iParticles) {
                Array<T, 3> pos(found[iParticle]->getPosition());
                pos = deltaX * pos + offset;
                Array<double, 3> doublePos(pos);
                fprintf(fp, "% .10e % .10e % .10e\n", doublePos[0], doublePos[1], doublePos[2]);
            }
        }
        fprintf(fp, "POINT_DATA %12lu\n", ulNumParticlesToWrite);
        std::map<plint, std::string>::const_iterator vectorIt = additionalVectors.begin();
        for (; vectorIt != additionalVectors.end(); ++vectorIt) {
            plint vectorID = vectorIt->first;
            std::string vectorName = vectorIt->second;
            fprintf(fp, "VECTORS ");
            fprintf(fp, "%s", vectorName.c_str());
            fprintf(fp, " double\n");
            if (numParticlesToWrite == numParticles) {
                for (pluint iParticle = 0; iParticle < numParticles; iParticle++) {
                    Array<T, 3> vectorValue;
                    found[iParticle]->getVector(vectorID, vectorValue);
                    Array<double, 3> doubleVectorValue(vectorValue);
                    fprintf(
                        fp, "% .10e % .10e % .10e\n", doubleVectorValue[0], doubleVectorValue[1],
                        doubleVectorValue[2]);
                }
            } else {
                for (pluint iParticle : iParticles) {
                    Array<T, 3> vectorValue;
                    found[iParticle]->getVector(vectorID, vectorValue);
                    Array<double, 3> doubleVectorValue(vectorValue);
                    fprintf(
                        fp, "% .10e % .10e % .10e\n", doubleVectorValue[0], doubleVectorValue[1],
                        doubleVectorValue[2]);
                }
            }
        }
        fprintf(fp, "SCALARS Tag double\n");
        fprintf(fp, "LOOKUP_TABLE default\n");
        if (numParticlesToWrite == numParticles) {
            for (pluint iParticle = 0; iParticle < numParticles; iParticle++) {
                double tag = (double)found[iParticle]->getTag();
                fprintf(fp, "% .10e\n", tag);
            }
        } else {
            for (pluint iParticle : iParticles) {
                double tag = (double)found[iParticle]->getTag();
                fprintf(fp, "% .10e\n", tag);
            }
        }
        std::map<plint, std::string>::const_iterator scalarIt = additionalScalars.begin();
        for (; scalarIt != additionalScalars.end(); ++scalarIt) {
            plint scalarID = scalarIt->first;
            std::string scalarName = scalarIt->second;
            fprintf(fp, "SCALARS ");
            fprintf(fp, "%s", scalarName.c_str());
            fprintf(fp, " double\n");
            fprintf(fp, "LOOKUP_TABLE default\n");
            if (numParticlesToWrite == numParticles) {
                for (pluint iParticle = 0; iParticle < numParticles; iParticle++) {
                    T scalarValue;
                    found[iParticle]->getScalar(scalarID, scalarValue);
                    double doubleScalarValue = (double)scalarValue;
                    fprintf(fp, "% .10e\n", doubleScalarValue);
                }
            } else {
                for (pluint iParticle : iParticles) {
                    T scalarValue;
                    found[iParticle]->getScalar(scalarID, scalarValue);
                    double doubleScalarValue = (double)scalarValue;
                    fprintf(fp, "% .10e\n", doubleScalarValue);
                }
            }
        }
        fclose(fp);
    }
}

template <typename T, template <typename U> class Descriptor>
void writeParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    T deltaX, Array<T, 3> const &offset, pluint maxNumParticlesToWrite)
{
    std::map<plint, std::string> additionalScalars;
    std::map<plint, std::string> additionalVectors;
    additionalVectors[0] = "Velocity";
    writeParticleVtk(
        particles, fName, additionalScalars, additionalVectors, deltaX, offset,
        maxNumParticlesToWrite);
}

template <typename T, template <typename U> class Descriptor>
void writeParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    std::map<plint, std::string> const &additionalScalars,
    std::map<plint, std::string> const &additionalVectors, T deltaX, Array<T, 3> const &offset,
    pluint maxNumParticlesToWrite)
{
    SparseBlockStructure3D blockStructure(particles.getBoundingBox());
    blockStructure.addBlock(particles.getBoundingBox(), 0);
    plint envelopeWidth = 1;
    MultiBlockManagement3D serialMultiBlockManagement(
        blockStructure, new OneToOneThreadAttribution, envelopeWidth);

    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > multiSerialParticles(
        serialMultiBlockManagement, defaultMultiBlockPolicy3D().getCombinedStatistics());

    copy(particles, particles.getBoundingBox(), multiSerialParticles, particles.getBoundingBox());
    particleVtkSerialImplementation(
        multiSerialParticles, fName, additionalScalars, additionalVectors, deltaX, offset,
        maxNumParticlesToWrite);
}

template <typename T, template <typename U> class Descriptor>
void writeSelectedParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    Box3D const &domain, util::SelectInt const &tags, T deltaX, Array<T, 3> const &offset)
{
    std::map<plint, std::string> additionalScalars;
    std::map<plint, std::string> additionalVectors;
    additionalVectors[0] = "Velocity";
    writeSelectedParticleVtk(
        particles, fName, domain, tags, additionalScalars, additionalVectors, deltaX, offset);
}

template <typename T, template <typename U> class Descriptor>
void writeSelectedParticleVtk(
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &particles, std::string const &fName,
    Box3D const &domain, util::SelectInt const &tags,
    std::map<plint, std::string> const &additionalScalars,
    std::map<plint, std::string> const &additionalVectors, T deltaX, Array<T, 3> const &offset)
{
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > filteredParticles(
        (MultiBlock3D &)particles);
    std::vector<MultiBlock3D *> particleArgs;
    particleArgs.push_back(&particles);
    particleArgs.push_back(&filteredParticles);
    applyProcessingFunctional(
        new CopySelectParticles3D<T, Descriptor>(tags.clone()), domain, particleArgs);
    SparseBlockStructure3D blockStructure(particles.getBoundingBox());
    blockStructure.addBlock(particles.getBoundingBox(), 0);
    plint envelopeWidth = 1;
    MultiBlockManagement3D serialMultiBlockManagement(
        blockStructure, new OneToOneThreadAttribution, envelopeWidth);

    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > multiSerialParticles(
        serialMultiBlockManagement, defaultMultiBlockPolicy3D().getCombinedStatistics());

    copy(filteredParticles, domain, multiSerialParticles, domain);
    particleVtkSerialImplementation(
        multiSerialParticles, fName, additionalScalars, additionalVectors, deltaX, offset, 0);
}

}  // namespace plb

#endif  // PARTICLE_VTK_3D_HH
