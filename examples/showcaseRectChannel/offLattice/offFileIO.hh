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

#ifndef OFF_FILE_IO_HH
#define OFF_FILE_IO_HH

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>

#include "core/runTimeDiagnostics.h"
#include "offLattice/connectedTriangleMesh.h"
#include "offLattice/connectedTriangleUtil.h"
#include "offLattice/offFileIO.h"
#include "offLattice/triangleMesh.h"

#define PLB_CBUFSIZ 4096  // Must be undefined at the end of this file.

namespace plb {

template <typename T>
OFFreader<T>::OFFreader(std::string fname)
{
    readOFF(fname);
}

template <typename T>
void OFFreader<T>::readOFF(std::string fname)
{
    FILE *fp = fopen(fname.c_str(), "rb");
    PLB_ASSERT(fp != 0);

    char buf[PLB_CBUFSIZ];
    char *sp = fgets(buf, PLB_CBUFSIZ, fp);
    if (sp == NULL) {
        throw PlbIOException("Problem reading file " + fname + ".");
    }

    char *cp = NULL;

    // Currently only ASCII files with header OFF can be read.
    cp = strstr(buf, "BINARY");
    if (cp != NULL) {
        throw PlbIOException("The file " + fname + " does not have an ASCII OFF formatting.");
    }
    std::vector<std::string> prefixes;
    prefixes.push_back("ST");
    prefixes.push_back("C");
    prefixes.push_back("N");
    prefixes.push_back("4");
    prefixes.push_back("n");
    for (int iPrefix = 0; iPrefix < (int)prefixes.size(); iPrefix++) {
        cp = strstr(buf, prefixes[iPrefix].c_str());
        if (cp != NULL) {
            throw PlbIOException("The file " + fname + " does not have a simple OFF formatting.");
        }
    }

    cp = strstr(buf, "OFF");

    bool failed = false;
    if (cp != NULL) {
        failed = readAsciiOFF(fp);
    } else {
        throw PlbIOException("The format of the file " + fname + " is not currently supported.");
    }
    fclose(fp);

    if (failed) {
        throw PlbIOException("Problem with file " + fname + ".");
    }
}

// Caution: in the readAsciiOFF method we keep all information that resides
//          in the OFF file. In the readAsciiOFF method of the TriangleSet class
//          we neglect all triangles that have one or more edges with length
//          equal to 0. We do not do this here, because from the OFFreader we
//          can directly create a connected mesh and we do not want to break
//          vertex and triangle numbering.
template <typename T>
bool OFFreader<T>::readAsciiOFF(FILE *fp)
{
    char commentCharacter = '#';
    long NVertices = 0, NFaces = 0, NEdges = 0;
    bool failed = false;
    if (readAhead(fp, commentCharacter) == EOF) {
        failed = true;
        return failed;
    }
    if (fscanf(fp, "%ld%ld%ld", &NVertices, &NFaces, &NEdges) != 3) {
        failed = true;
        return failed;
    }

    char fmt[32];
    if (sizeof(T) == sizeof(float)) {
        strcpy(fmt, "%f%f%f");
    } else if (sizeof(T) == sizeof(double)) {
        strcpy(fmt, "%lf%lf%lf");
    } else if (sizeof(T) == sizeof(long double)) {
        strcpy(fmt, "%Lf%Lf%Lf");
    } else {
        failed = true;
        return failed;
    }

    vertices.clear();
    vertices.resize((size_t)NVertices);
    for (long iVertex = 0; iVertex < NVertices; iVertex++) {
        if (readAhead(fp, commentCharacter) == EOF) {
            failed = true;
            return failed;
        }
        if (fscanf(fp, fmt, &vertices[iVertex][0], &vertices[iVertex][1], &vertices[iVertex][2])
            != 3) {
            failed = true;
            return failed;
        }
    }

    facets.clear();
    facets.resize((size_t)NFaces);
    for (long iFace = 0; iFace < NFaces; iFace++) {
        if (readAhead(fp, commentCharacter) == EOF) {
            failed = true;
            return failed;
        }
        long Nv;
        if (fscanf(fp, "%ld", &Nv) != 1) {
            failed = true;
            return failed;
        }
        facets[iFace].resize((size_t)Nv);
        for (long iVertex = 0; iVertex < Nv; iVertex++) {
            long ind;
            if (fscanf(fp, "%ld", &ind) != 1) {
                failed = true;
                return failed;
            }
            if (ind < 0 || ind >= NVertices) {
                failed = true;
                return failed;
            }
            facets[iFace][iVertex] = ind;
        }
    }

    return failed;
}

/// Skip nLines number of lines in the file.
template <typename T>
void OFFreader<T>::skipLines(plint nLines, FILE *fp) const
{
    for (plint i = 0; i < nLines; i++)
        while (fgetc(fp) != '\n')
            ;
}

/// Read the file character-by-character. This function returns 0 if a non-white-space
/// character is found, or EOF if the end-of-file is found. It "consumes" the
/// rest of the line if the "commentCharacter" is found. Obviously, it works only
/// with text files.
template <typename T>
int OFFreader<T>::readAhead(FILE *fp, char commentCharacter) const
{
    int nextChar;
    while ((nextChar = fgetc(fp)) != EOF) {
        if (isspace(nextChar)) {
            continue;
        } else if (nextChar == commentCharacter) {
            skipLines(1, fp);
        } else {
#ifdef PLB_DEBUG
            int rv = ungetc(nextChar, fp);
#endif
            PLB_ASSERT(rv != EOF);  // Unexpected error.
            return 0;
        }
    }
    return EOF;
}

template <typename T>
void writeAsciiOFF(TriangleMesh<T> &mesh, std::string fname, T eps, int numDecimalDigits)
{
    typedef typename ConnectedTriangleMesh<T>::PTriangleIterator PTriangleIterator;
    typedef typename ConnectedTriangleMesh<T>::PVertexIterator PVertexIterator;
    typedef typename ConnectedTriangleMesh<T>::PTriangle PTriangle;
    typedef typename ConnectedTriangleMesh<T>::PVertex PVertex;

    RawConnectedTriangleMesh<T> connectedMesh = MeshConnector<T>(mesh, eps).generateConnectedMesh();

    if (global::mpi().isMainProcessor()) {
        if (connectedMesh.getNumTriangles() == 0) {
            return;
        }
        std::ofstream ofile(fname.c_str());
        if (ofile.is_open()) {
            ofile << "OFF\n";
            ofile << connectedMesh.getNumVertices() << " " << connectedMesh.getNumTriangles() << " "
                  << 0 << '\n';

            plint vertexIDtagging = connectedMesh.getVertexTag("UniqueID");
            PVertexIterator vertexIterator(connectedMesh.vertexIterator());
            std::map<plint, plint> toNewVertexID;
            plint newVertexID = 0;
            while (!vertexIterator->end()) {
                PVertex vertex(vertexIterator->next());
                Array<T, 3> v = vertex->get();
                ofile << std::setprecision(numDecimalDigits) << std::scientific << v[0] << " "
                      << v[1] << " " << v[2] << '\n';
                toNewVertexID[vertex->tag(vertexIDtagging)] = newVertexID;
                ++newVertexID;
            }

            PTriangleIterator triangleIterator(connectedMesh.triangleIterator());
            while (!triangleIterator->end()) {
                PTriangle triangle(triangleIterator->next());
                PVertex v0, v1, v2;
                v0 = triangle->vertex(0);
                v1 = triangle->vertex(1);
                v2 = triangle->vertex(2);
                plint ind0, ind1, ind2;
                ind0 = v0->tag(vertexIDtagging);
                ind1 = v1->tag(vertexIDtagging);
                ind2 = v2->tag(vertexIDtagging);
                ofile << 3 << " " << toNewVertexID[ind0] << " " << toNewVertexID[ind1] << " "
                      << toNewVertexID[ind2] << '\n';
            }
        }
    }
}

}  // namespace plb

#undef PLB_CBUFSIZ

#endif  // OFF_FILE_IO_HH
