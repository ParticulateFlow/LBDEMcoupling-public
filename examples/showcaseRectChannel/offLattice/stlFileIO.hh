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

#ifndef STL_FILE_IO_HH
#define STL_FILE_IO_HH

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include "core/runTimeDiagnostics.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/stlFileIO.h"
#include "offLattice/triangleMesh.h"

#define PLB_CBUFSIZ \
    4096  // Must be undefined at the end of this file.
          // Must be larger than 80 which is the size of the binary STL header.
          // Must be large enough to contain all characters of an ASCII STL
          // file including at least the first "endfacet" word.

namespace plb {

template <typename T>
STLreader<T>::STLreader(std::string fname, T eps_, TriangleSelector<T> *selector) : eps(eps_)
{
    readSTL(fname, selector);
    delete selector;
}

template <typename T>
void STLreader<T>::readSTL(std::string fname, TriangleSelector<T> *selector)
{
    FILE *fp = fopen(fname.c_str(), "rb");
    PLB_ASSERT(fp != 0);

    bool failed = false;
    if (isAsciiSTL(fp, fname)) {
        failed = readAsciiSTL(fp, selector);
    } else {
        failed = readBinarySTL(fp, selector);
    }
    fclose(fp);

    if (failed) {
        throw PlbIOException("Problem with file " + fname + ".");
    }
}

template <typename T>
bool STLreader<T>::isAsciiSTL(FILE *fp, std::string fname)
{
    char buf[PLB_CBUFSIZ + 1];

#ifdef PLB_DEBUG
    size_t sz = fread(buf, sizeof(char), PLB_CBUFSIZ, fp);
    if (sz != PLB_CBUFSIZ) {
        PLB_ASSERT(!ferror(fp));
    }
#else
    (void)fread(buf, sizeof(char), PLB_CBUFSIZ, fp);
#endif
    buf[PLB_CBUFSIZ] = '\0';
    rewind(fp);

    if (strstr(buf, "solid") == NULL) {  // If "solid" does not exist then it is binary STL.
        return false;
    }

    int rv = fseek(fp, 80L, SEEK_SET);
    if (rv == -1) {
        throw PlbIOException("Problem seeking the file " + fname + ".");
    }

#ifdef PLB_DEBUG
    sz = fread(buf, sizeof(char), PLB_CBUFSIZ, fp);
    if (sz != PLB_CBUFSIZ) {
        PLB_ASSERT(!ferror(fp));
    }
#else
    (void)fread(buf, sizeof(char), PLB_CBUFSIZ, fp);
#endif
    buf[PLB_CBUFSIZ] = '\0';
    rewind(fp);

    if (strstr(buf, "endfacet") != NULL) {  // If "endfacet" exists then it is ASCII STL.
        return true;
    }

    return false;
}

template <typename T>
bool STLreader<T>::readAsciiSTL(FILE *fp, TriangleSelector<T> *selector)
{
    char buf[PLB_CBUFSIZ];
    char *cp, *sp;
    bool failed = false;

    sp = fgets(buf, PLB_CBUFSIZ, fp);
    if (sp == NULL) {
        failed = true;
        return failed;
    }
    if (!checkForBufferOverflow(buf)) {
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

    parts.clear();
    names.clear();

    cp = strstr(buf, "solid");
    if (cp == NULL) {
        failed = true;
        return failed;
    }
    char *namep = cp + 5;
    namep = cleanString(namep);
    names.push_back(std::string(namep));

    while (cp != NULL) {
        if (fgets(buf, PLB_CBUFSIZ, fp) == NULL) {
            failed = true;
            return failed;
        }
        if (!checkForBufferOverflow(buf)) {
            failed = true;
            return failed;
        }

        if ((cp = strstr(buf, "color")) != NULL) {
            if (fgets(buf, PLB_CBUFSIZ, fp) == NULL) {
                failed = true;
                return failed;
            }
            if (!checkForBufferOverflow(buf)) {
                failed = true;
                return failed;
            }
        }

        std::vector<Triangle> triangles;
        do {
            if ((cp = strstr(buf, "facet normal")) == NULL) {
                failed = true;
                return failed;
            }
            cp += 12;
            Array<T, 3> n;
            if (sscanf(cp, fmt, &n[0], &n[1], &n[2]) != 3) {
                failed = true;
                return failed;
            }

            if (fgets(buf, PLB_CBUFSIZ, fp) == NULL || strstr(buf, "outer loop") == NULL) {
                failed = true;
                return failed;
            }
            if (!checkForBufferOverflow(buf)) {
                failed = true;
                return failed;
            }

            Triangle triangle;
            for (int i = 0; i < 3; i++) {
                if (fgets(buf, PLB_CBUFSIZ, fp) == NULL || (cp = strstr(buf, "vertex")) == NULL) {
                    failed = true;
                    return failed;
                }
                if (!checkForBufferOverflow(buf)) {
                    failed = true;
                    return failed;
                }
                cp += 6;
                triangle[i][0] = T();
                triangle[i][1] = T();
                triangle[i][2] = T();
                if (sscanf(cp, fmt, &triangle[i][0], &triangle[i][1], &triangle[i][2]) != 3) {
                    failed = true;
                    return failed;
                }
            }

            if (fgets(buf, PLB_CBUFSIZ, fp) == NULL || strstr(buf, "endloop") == NULL) {
                failed = true;
                return failed;
            }
            if (!checkForBufferOverflow(buf)) {
                failed = true;
                return failed;
            }

            if (fgets(buf, PLB_CBUFSIZ, fp) == NULL || strstr(buf, "endfacet") == NULL) {
                failed = true;
                return failed;
            }
            if (!checkForBufferOverflow(buf)) {
                failed = true;
                return failed;
            }

            // We keep triangles that have a zero area, since they might be used to
            // fix topological problems with the connectivity of the triangle set.
            // However, we do not keep triangles that have one or more edges with a
            // zero length. This has as a result, that if a triangle is kept in the
            // triangle set and has zero area, its vertices are distinct but belong
            // to a straight line.
            if (!triangleHasZeroLengthEdges(triangle)) {
                fixOrientation(triangle, n);
                plint partId = (plint)parts.size();
                if (selector == 0 || (*selector)(triangle, partId)) {
                    triangles.push_back(triangle);
                }
            }

            if (fgets(buf, PLB_CBUFSIZ, fp) == NULL) {
                failed = true;
                return failed;
            }
            cp = strstr(buf, "endsolid");
            if (cp == NULL) {
                if (!checkForBufferOverflow(buf)) {
                    failed = true;
                    return failed;
                }
            }
        } while (cp == NULL);

        parts.push_back(triangles);

        if (fgets(buf, PLB_CBUFSIZ, fp) == NULL)
            break;

        cp = strstr(buf, "solid");

        if (cp != NULL) {
            if (!checkForBufferOverflow(buf)) {
                failed = true;
                return failed;
            }
            namep = cp + 5;
            namep = cleanString(namep);
            names.push_back(std::string(namep));
        }
    }

    if (parts.size() != names.size()) {
        failed = true;
    }

    return failed;
}

template <typename T>
bool STLreader<T>::readBinarySTL(FILE *fp, TriangleSelector<T> *selector)
{
    char buf[PLB_CBUFSIZ];
    std::fill(buf, buf + PLB_CBUFSIZ, '\0');
    unsigned int nt;
    float array[3];
    unsigned short abc;
    bool failed = false;

    parts.clear();
    names.clear();

    int count = 0;
    while (fread(buf, sizeof(char), 80, fp) == 80 && fread(&nt, sizeof(unsigned int), 1, fp) == 1) {
        count++;

        buf[81] = '\0';
        char *namep = cleanString(buf);
        names.push_back(std::string(namep));

        std::vector<Triangle> triangles;
        for (unsigned int it = 0; it < nt; it++) {
            if (fread(array, sizeof(float), 3, fp) != 3) {
                failed = true;
                return failed;
            }
            Array<T, 3> n;
            n[0] = array[0];
            n[1] = array[1];
            n[2] = array[2];

            Triangle triangle;
            for (int i = 0; i < 3; i++) {
                if (fread(array, sizeof(float), 3, fp) != 3) {
                    failed = true;
                    return failed;
                }
                triangle[i][0] = array[0];
                triangle[i][1] = array[1];
                triangle[i][2] = array[2];
            }

            if (fread(&abc, sizeof(unsigned short), 1, fp) != 1) {
                failed = true;
                return failed;
            }

            // We keep triangles that have a zero area, since they might be used to
            // fix topological problems with the connectivity of the triangle set.
            // However, we do not keep triangles that have one or more edges with a
            // zero length. This has as a result, that if a triangle is kept in the
            // triangle set and has zero area, its vertices are distinct but belong
            // to a straight line.
            if (!triangleHasZeroLengthEdges(triangle)) {
                fixOrientation(triangle, n);
                plint partId = (plint)parts.size();
                if (selector == 0 || (*selector)(triangle, partId)) {
                    triangles.push_back(triangle);
                }
            }
        }
        parts.push_back(triangles);
    }

    if (count == 0 || parts.size() != names.size()) {
        failed = true;
        return failed;
    }

    return failed;
}

template <typename T>
bool STLreader<T>::triangleHasZeroLengthEdges(typename STLreader<T>::Triangle const &triangle) const
{
    if (util::isZero(norm(triangle[1] - triangle[0]), eps)
        || util::isZero(norm(triangle[2] - triangle[0]), eps)
        || util::isZero(norm(triangle[2] - triangle[1]), eps))
    {
        return true;
    }
    return false;
}

template <typename T>
void STLreader<T>::fixOrientation(
    typename STLreader<T>::Triangle &triangle, Array<T, 3> const &n) const
{
    bool isAreaWeighted = false;
    Array<T, 3> computedNormal =
        computeTriangleNormal(triangle[0], triangle[1], triangle[2], isAreaWeighted);
    if (dot(computedNormal, n) < (T)0) {
        std::swap(triangle[1], triangle[2]);
    }
}

/// Check to see if the buffer that holds a line of text is full or not.
template <typename T>
bool STLreader<T>::checkForBufferOverflow(char *buf)
{
    while (*buf != '\0') {
        if (*buf == '\n') {
            return true;  // All the line was read into the buffer.
        }
        buf++;
    }
    return false;  // The full line of text was not read into the buffer.
}

/// Clean a null-terminated string from leading and trailing white space.
template <typename T>
char *STLreader<T>::cleanString(char *str)
{
    size_t len = strlen(str);
    if (len == 0) {
        return str;
    }

    char *cp = str + len - 1;
    while ((cp + 1) != str && isspace(*cp--)) {
        *(cp + 1) = '\0';
    }

    cp = str;
    while (*cp && isspace(*cp++)) {
        ;
    }
    if (cp != str) {
        cp--;
    }

    return cp;
}

template <typename T>
void writeAsciiSTL(
    TriangleMesh<T> &mesh, std::string fname, int numDecimalDigits, bool mainProcOnly)
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        {
            typename TriangleMesh<T>::PTriangleIterator iterator(mesh.triangleIterator());
            plint numTriangles = 0;
            while (!iterator->end()) {
                iterator->next();
                numTriangles++;
            }
            if (numTriangles == 0) {
                return;
            }
        }

        typename TriangleMesh<T>::PTriangleIterator iterator(mesh.triangleIterator());
        std::ofstream ofile(fname.c_str());
        if (ofile.is_open()) {
            ofile << "solid plb\n";
            while (!iterator->end()) {
                typename TriangleMesh<T>::PTriangle triangle(iterator->next());
                Array<T, 3> v0 = (*triangle)[0];
                Array<T, 3> v1 = (*triangle)[1];
                Array<T, 3> v2 = (*triangle)[2];

                Array<T, 3> n = triangle->normal();
                ofile << "  facet normal " << std::setprecision(numDecimalDigits) << std::scientific
                      << n[0] << " " << n[1] << " " << n[2] << "\n";
                ofile << "    outer loop\n";
                ofile << "      vertex " << std::setprecision(numDecimalDigits) << std::scientific
                      << v0[0] << " " << v0[1] << " " << v0[2] << "\n";
                ofile << "      vertex " << std::setprecision(numDecimalDigits) << std::scientific
                      << v1[0] << " " << v1[1] << " " << v1[2] << "\n";
                ofile << "      vertex " << std::setprecision(numDecimalDigits) << std::scientific
                      << v2[0] << " " << v2[1] << " " << v2[2] << "\n";
                ofile << "    endloop\n";
                ofile << "  endfacet\n";
            }
            ofile << "endsolid plb\n";
        }
    }
}

template <typename T>
void writeMultiPartAsciiSTL(
    TriangleMesh<T> &mesh, std::string fname, int numDecimalDigits, bool mainProcOnly)
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        {
            typename TriangleMesh<T>::PTriangleIterator iterator(mesh.triangleIterator());
            plint numTriangles = 0;
            while (!iterator->end()) {
                iterator->next();
                numTriangles++;
            }
            if (numTriangles == 0) {
                return;
            }
        }

        std::ofstream ofile(fname.c_str());
        if (ofile.is_open()) {
            for (plint iPart = 0; iPart < mesh.numParts(); iPart++) {
                typename TriangleMesh<T>::PTriangleIterator iterator(mesh.triangleIterator(iPart));
                ofile << "solid " << mesh.partName(iPart) << "\n";
                while (!iterator->end()) {
                    typename TriangleMesh<T>::PTriangle triangle(iterator->next());
                    Array<T, 3> v0 = (*triangle)[0];
                    Array<T, 3> v1 = (*triangle)[1];
                    Array<T, 3> v2 = (*triangle)[2];

                    Array<T, 3> n = triangle->normal();
                    ofile << "  facet normal " << std::setprecision(numDecimalDigits)
                          << std::scientific << n[0] << " " << n[1] << " " << n[2] << "\n";
                    ofile << "    outer loop\n";
                    ofile << "      vertex " << std::setprecision(numDecimalDigits)
                          << std::scientific << v0[0] << " " << v0[1] << " " << v0[2] << "\n";
                    ofile << "      vertex " << std::setprecision(numDecimalDigits)
                          << std::scientific << v1[0] << " " << v1[1] << " " << v1[2] << "\n";
                    ofile << "      vertex " << std::setprecision(numDecimalDigits)
                          << std::scientific << v2[0] << " " << v2[1] << " " << v2[2] << "\n";
                    ofile << "    endloop\n";
                    ofile << "  endfacet\n";
                }
                ofile << "endsolid " << mesh.partName(iPart) << "\n";
            }
        }
    }
}

template <typename T>
void writeBinarySTL(
    typename TriangleMesh<T>::PTriangleIterator const &iterator, std::string fname,
    bool mainProcOnly)
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        typename TriangleMesh<T>::PTriangleIterator it2(iterator->clone());
        plint numTriangles = 0;
        while (!it2->end()) {
            it2->next();
            numTriangles++;
        }
        if (numTriangles == 0) {
            return;
        }

        FILE *fp = fopen(fname.c_str(), "wb");
        if (fp == 0) {
            return;
        }

        unsigned int nt = numTriangles;
        unsigned short abc = 0;
        char buf[80] = {'\0'};

        strcpy(buf, "plb");

        fwrite(buf, sizeof(char), 80, fp);
        fwrite(&nt, sizeof(unsigned int), 1, fp);
        for (unsigned int i = 0; i < nt; i++) {
            typename TriangleMesh<T>::PTriangle triangle(iterator->next());

            Array<T, 3> v0 = (*triangle)[0];
            Array<T, 3> v1 = (*triangle)[1];
            Array<T, 3> v2 = (*triangle)[2];

            Array<T, 3> nrml = triangle->normal();

            float n[3];
            n[0] = nrml[0];
            n[1] = nrml[1];
            n[2] = nrml[2];
            fwrite((void *)n, sizeof(float), 3, fp);
            float v[3];
            v[0] = v0[0];
            v[1] = v0[1];
            v[2] = v0[2];
            fwrite((void *)v, sizeof(float), 3, fp);
            v[0] = v1[0];
            v[1] = v1[1];
            v[2] = v1[2];
            fwrite((void *)v, sizeof(float), 3, fp);
            v[0] = v2[0];
            v[1] = v2[1];
            v[2] = v2[2];
            fwrite((void *)v, sizeof(float), 3, fp);
            fwrite(&abc, sizeof(unsigned short), 1, fp);
        }

        fclose(fp);
    }
}

template <typename T>
void writeBinarySTL(TriangleMesh<T> &mesh, std::string fname, bool mainProcOnly)
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        writeBinarySTL<T>(mesh.triangleIterator(), fname, mainProcOnly);
    }
}

template <typename T>
void writeMultiPartBinarySTL(TriangleMesh<T> &mesh, std::string fname, bool mainProcOnly)
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        {
            typename TriangleMesh<T>::PTriangleIterator iterator(mesh.triangleIterator());
            plint numTriangles = 0;
            while (!iterator->end()) {
                iterator->next();
                numTriangles++;
            }
            if (numTriangles == 0) {
                return;
            }
        }
        FILE *fp = fopen(fname.c_str(), "wb");
        if (fp == 0) {
            return;
        }

        for (plint iPart = 0; iPart < mesh.numParts(); iPart++) {
            typename TriangleMesh<T>::PTriangleIterator iterator(mesh.triangleIterator(iPart));

            plint numTriangles = 0;
            while (!iterator->end()) {
                iterator->next();
                numTriangles++;
            }

            unsigned int nt = numTriangles;
            unsigned short abc = 0;
            char buf[80] = {'\0'};

            char *name = (char *)malloc((mesh.partName(iPart).length() + 1) * sizeof(char));
            strcpy(name, mesh.partName(iPart).c_str());
            if (strlen(name) > 79) {
                name[79] = '\0';
            }
            strcpy(buf, name);
            free(name);

            fwrite(buf, sizeof(char), 80, fp);
            fwrite(&nt, sizeof(unsigned int), 1, fp);
            iterator = mesh.triangleIterator(iPart);
            for (unsigned int i = 0; i < nt; i++) {
                typename TriangleMesh<T>::PTriangle triangle(iterator->next());

                Array<T, 3> v0 = (*triangle)[0];
                Array<T, 3> v1 = (*triangle)[1];
                Array<T, 3> v2 = (*triangle)[2];

                Array<T, 3> nrml = triangle->normal();

                float n[3];
                n[0] = nrml[0];
                n[1] = nrml[1];
                n[2] = nrml[2];
                fwrite((void *)n, sizeof(float), 3, fp);
                float v[3];
                v[0] = v0[0];
                v[1] = v0[1];
                v[2] = v0[2];
                fwrite((void *)v, sizeof(float), 3, fp);
                v[0] = v1[0];
                v[1] = v1[1];
                v[2] = v1[2];
                fwrite((void *)v, sizeof(float), 3, fp);
                v[0] = v2[0];
                v[1] = v2[1];
                v[2] = v2[2];
                fwrite((void *)v, sizeof(float), 3, fp);
                fwrite(&abc, sizeof(unsigned short), 1, fp);
            }
        }

        fclose(fp);
    }
}

}  // namespace plb

#undef PLB_CBUFSIZ

#endif  // STL_FILE_IO_HH
