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

#ifndef VTK_FILE_IO_HH
#define VTK_FILE_IO_HH

#include <fstream>
#include <limits>

#include "io/imageWriter.h"
#include "offLattice/connectedTriangleUtil.h"
#include "offLattice/vtkFileIO.h"

namespace plb {

template <typename T>
void writeMultiPartVTK(
    RawConnectedTriangleMesh<T> &mesh, FileName fname, bool writeVertexNormals,
    std::string vertexNormalsName)
{
    if (!global::mpi().isMainProcessor())
        return;
    if (mesh.getNumTriangles() == 0) {
        pcout << "Warning: Trying to save a mesh with no triangles." << std::endl;
        return;
    }

    for (plint iPart = 0; iPart < mesh.numParts(); ++iPart) {
        FileName partFileName(fname);
        if (mesh.numParts() > 1) {
            partFileName.setName(createFileName(partFileName.getName() + "_", iPart, 3));
        }
        RawConnectedTriangleMesh<T> part = extractConnectedPart(mesh, iPart);
        writeVTK(part, partFileName, writeVertexNormals, vertexNormalsName);
    }
}

template <typename T>
void writeVTK(
    RawConnectedTriangleMesh<T> &mesh, FileName fname, bool writeVertexNormals,
    std::string vertexNormalsName)
{
    if (!global::mpi().isMainProcessor())
        return;
    if (mesh.getNumTriangles() == 0)
        return;
    typedef typename ConnectedTriangleMesh<T>::PTriangleIterator PTriangleIterator;
    typedef typename ConnectedTriangleMesh<T>::PVertexIterator PVertexIterator;
    typedef typename ConnectedTriangleMesh<T>::PTriangle PTriangle;
    typedef typename ConnectedTriangleMesh<T>::PVertex PVertex;

    std::ofstream ofile(fname.get().c_str());
    ofile.precision(10);
    std::scientific(ofile);

    ofile << "# vtk DataFile Version 3.0\n";
    ofile << "Surface mesh created with Palabos\n";
    ofile << "ASCII\n";
    ofile << "DATASET UNSTRUCTURED_GRID\n";

    ofile << "POINTS " << mesh.getNumVertices()
          << (sizeof(T) == sizeof(double) ? " double" : " float") << "\n";

    plint vertexIDtagging = mesh.getVertexTag("UniqueID");
    PVertexIterator vertexIt = mesh.vertexIterator();
    std::map<plint, plint> toNewVertexID;
    plint newVertexID = 0;
    while (!vertexIt->end()) {
        PVertex vertex = vertexIt->next();
        ofile << (*vertex)[0] << " " << (*vertex)[1] << " " << (*vertex)[2] << "\n";
        toNewVertexID[vertex->tag(vertexIDtagging)] = newVertexID;
        ++newVertexID;
    }
    ofile << "\n";

    ofile << "CELLS " << mesh.getNumTriangles() << " " << 4 * mesh.getNumTriangles() << "\n";

    PTriangleIterator triangleIt = mesh.triangleIterator();
    while (!triangleIt->end()) {
        PTriangle triangle = triangleIt->next();
        plint i0 = triangle->vertex(0)->tag(vertexIDtagging);
        plint i1 = triangle->vertex(1)->tag(vertexIDtagging);
        plint i2 = triangle->vertex(2)->tag(vertexIDtagging);
        ofile << "3 " << i0 << " " << i1 << " " << i2 << "\n";
    }
    ofile << "\n";

    ofile << "CELL_TYPES " << mesh.getNumTriangles() << "\n";

    for (plint i = 0; i < mesh.getNumTriangles(); ++i) {
        ofile << "5\n";
    }
    ofile << "\n";

    ofile << "POINT_DATA " << mesh.getNumVertices() << "\n";

    for (plint property = 0; property < mesh.numVertexProperties(); ++property) {
        ofile << "SCALARS " << mesh.getVertexPropertyName(property)
              << (sizeof(T) == sizeof(double) ? " double" : " float") << " 1\n"
              << "LOOKUP_TABLE default\n";
        vertexIt = mesh.vertexIterator();
        while (!vertexIt->end()) {
            PVertex vertex = vertexIt->next();
            ofile << vertex->property(property) << "\n";
        }
        ofile << "\n";
    }

    if (writeVertexNormals) {
        ofile << "VECTORS " << vertexNormalsName
              << (sizeof(T) == sizeof(double) ? " double" : " float") << "\n";
        vertexIt = mesh.vertexIterator();
        while (!vertexIt->end()) {
            PVertex vertex = vertexIt->next();
            Array<T, 3> n = vertex->normal();
            ofile << n[0] << " " << n[1] << " " << n[2] << "\n";
        }
        ofile << "\n";
    }
}

}  // namespace plb

#endif  // VTK_FILE_IO_HH
