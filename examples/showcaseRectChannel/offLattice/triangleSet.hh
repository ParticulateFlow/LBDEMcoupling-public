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

#ifndef TRIANGLE_SET_HH
#define TRIANGLE_SET_HH

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "offLattice/triangleSet.h"

#define PLB_CBUFSIZ \
    4096  // Must be undefined at the end of this file.
          // Must be larger than 80 which is the size of the binary STL header.
          // Must be large enough to contain all characters of an ASCII STL
          // file including at least the first "endfacet" word.

namespace plb {

template <typename T>
TriangleSet<T>::TriangleSet(Precision precision_) :
    minEdgeLength(std::numeric_limits<T>::max()),
    maxEdgeLength(std::numeric_limits<T>::min()),
    minTriangleArea(std::numeric_limits<T>::max()),
    maxTriangleArea(std::numeric_limits<T>::min())
{
    PLB_ASSERT(precision_ == FLT || precision_ == DBL || precision_ == LDBL || precision_ == INF);
    eps = plb::getEpsilon<T>(precision_);

    T huge = std::numeric_limits<T>::max();

    boundingCuboid.lowerLeftCorner = Array<T, 3>((T)huge, (T)huge, (T)huge);
    boundingCuboid.upperRightCorner = Array<T, 3>((T)-huge, (T)-huge, (T)-huge);
}

template <typename T>
TriangleSet<T>::TriangleSet(T eps_) :
    minEdgeLength(std::numeric_limits<T>::max()),
    maxEdgeLength(std::numeric_limits<T>::min()),
    minTriangleArea(std::numeric_limits<T>::max()),
    maxTriangleArea(std::numeric_limits<T>::min())
{
    PLB_ASSERT(eps_ > (T)0);
    eps = eps_;

    T huge = std::numeric_limits<T>::max();

    boundingCuboid.lowerLeftCorner = Array<T, 3>((T)huge, (T)huge, (T)huge);
    boundingCuboid.upperRightCorner = Array<T, 3>((T)-huge, (T)-huge, (T)-huge);
}

template <typename T>
TriangleSet<T>::TriangleSet(
    std::vector<Triangle> const &triangles_, Precision precision_, TriangleSelector<T> *selector) :
    minEdgeLength(std::numeric_limits<T>::max()),
    maxEdgeLength(std::numeric_limits<T>::min()),
    minTriangleArea(std::numeric_limits<T>::max()),
    maxTriangleArea(std::numeric_limits<T>::min())
{
    PLB_ASSERT(precision_ == FLT || precision_ == DBL || precision_ == LDBL || precision_ == INF);
    eps = plb::getEpsilon<T>(precision_);

    plint partId = 0;  // There is only one part.

    for (pluint i = 0; i < triangles_.size(); ++i) {
        // We keep triangles that have a zero area, since they might be used to
        // fix topological problems with the connectivity of the triangle set.
        // However, we do not keep triangles that have one or more edges with a
        // zero length. This has as a result, that if a triangle is kept in the
        // triangle set and has zero area, its vertices are distinct but belong
        // to a straight line.
        if (!triangleHasZeroLengthEdges(triangles_[i], eps)) {
            if (selector == 0 || (*selector)(triangles_[i], partId)) {
                triangles.push_back(triangles_[i]);
            }
        }
    }

    delete selector;

    computeMinMaxEdges();
    computeMinMaxAreas();
    computeBoundingCuboid();
}

template <typename T>
TriangleSet<T>::TriangleSet(
    std::vector<Triangle> const &triangles_, T eps_, TriangleSelector<T> *selector) :
    minEdgeLength(std::numeric_limits<T>::max()),
    maxEdgeLength(std::numeric_limits<T>::min()),
    minTriangleArea(std::numeric_limits<T>::max()),
    maxTriangleArea(std::numeric_limits<T>::min())
{
    PLB_ASSERT(eps_ > (T)0);
    eps = eps_;

    plint partId = 0;  // There is only one part.

    for (pluint i = 0; i < triangles_.size(); ++i) {
        // We keep triangles that have a zero area, since they might be used to
        // fix topological problems with the connectivity of the triangle set.
        // However, we do not keep triangles that have one or more edges with a
        // zero length. This has as a result, that if a triangle is kept in the
        // triangle set and has zero area, its vertices are distinct but belong
        // to a straight line.
        if (!triangleHasZeroLengthEdges(triangles_[i], eps)) {
            if (selector == 0 || (*selector)(triangles_[i], partId)) {
                triangles.push_back(triangles_[i]);
            }
        }
    }

    delete selector;

    computeMinMaxEdges();
    computeMinMaxAreas();
    computeBoundingCuboid();
}

template <typename T>
TriangleSet<T>::TriangleSet(
    std::string fname, Precision precision_, SurfaceGeometryFileFormat fformat,
    TriangleSelector<T> *selector) :
    minEdgeLength(std::numeric_limits<T>::max()),
    maxEdgeLength(std::numeric_limits<T>::min()),
    minTriangleArea(std::numeric_limits<T>::max()),
    maxTriangleArea(std::numeric_limits<T>::min())
{
    PLB_ASSERT(precision_ == FLT || precision_ == DBL || precision_ == LDBL || precision_ == INF);
    PLB_ASSERT(fformat == STL || fformat == OFF);

    eps = plb::getEpsilon<T>(precision_);

    switch (fformat) {
    case STL:
    default:
        readSTL(fname, selector);
        break;
    case OFF:
        readOFF(fname, selector);
        break;
    }

    delete selector;

    computeMinMaxEdges();
    computeMinMaxAreas();
    computeBoundingCuboid();
}

template <typename T>
TriangleSet<T>::TriangleSet(
    std::string fname, T eps_, SurfaceGeometryFileFormat fformat, TriangleSelector<T> *selector) :
    minEdgeLength(std::numeric_limits<T>::max()),
    maxEdgeLength(std::numeric_limits<T>::min()),
    minTriangleArea(std::numeric_limits<T>::max()),
    maxTriangleArea(std::numeric_limits<T>::min())
{
    PLB_ASSERT(eps_ > (T)0);
    PLB_ASSERT(fformat == STL || fformat == OFF);

    eps = eps_;

    switch (fformat) {
    case STL:
    default:
        readSTL(fname, selector);
        break;
    case OFF:
        readOFF(fname, selector);
        break;
    }

    delete selector;

    computeMinMaxEdges();
    computeMinMaxAreas();
    computeBoundingCuboid();
}

template <typename T>
std::vector<typename TriangleSet<T>::Triangle> const &TriangleSet<T>::getTriangles() const
{
    return triangles;
}

template <typename T>
void TriangleSet<T>::setPrecision(Precision precision_)
{
    PLB_ASSERT(precision_ == FLT || precision_ == DBL || precision_ == LDBL || precision_ == INF);
    eps = plb::getEpsilon<T>(precision_);
}

template <typename T>
void TriangleSet<T>::setEpsilon(T eps_)
{
    PLB_ASSERT(eps_ > (T)0);
    eps = eps_;
}

template <typename T>
void TriangleSet<T>::readSTL(std::string fname, TriangleSelector<T> *selector)
{
    FILE *fp = fopen(fname.c_str(), "rb");
    PLB_ASSERT(fp != 0);  // The input file cannot be read.

    if (isAsciiSTL(fp)) {
        readAsciiSTL(fp, selector);
    } else {
        readBinarySTL(fp, selector);
    }

    fclose(fp);
}

template <typename T>
bool TriangleSet<T>::isAsciiSTL(FILE *fp)
{
    char buf[PLB_CBUFSIZ + 1];

#ifdef PLB_DEBUG
    size_t sz = fread(buf, sizeof(char), PLB_CBUFSIZ, fp);
    if (sz != PLB_CBUFSIZ) {
        PLB_ASSERT(!ferror(fp));
    }
#else
    (void)fread(buf, sizeof(char), PLB_CBUFSIZ, fp);  // TODO: (void) does not suppress warning. We
                                                      // should alway check return. throw exception.
#endif
    buf[PLB_CBUFSIZ] = '\0';
    rewind(fp);

    if (strstr(buf, "solid") == NULL) {  // If "solid" does not exist then it is binary STL.
        return false;
    }

#ifdef PLB_DEBUG
    int rv = fseek(fp, 80L, SEEK_SET);
#else
    (void)fseek(fp, 80L, SEEK_SET);
#endif
    PLB_ASSERT(rv != -1);

#ifdef PLB_DEBUG
    sz = fread(buf, sizeof(char), PLB_CBUFSIZ, fp);
    if (sz != PLB_CBUFSIZ) {
        PLB_ASSERT(!ferror(fp));
    }
#else
    (void)fread(buf, sizeof(char), PLB_CBUFSIZ, fp);  // TODO: (void) does not suppress warning. We
                                                      // should alway check return. throw exception.
#endif
    buf[PLB_CBUFSIZ] = '\0';
    rewind(fp);

    if (strstr(buf, "endfacet") != NULL) {  // If "endfacet" exists then it is ASCII STL.
        return true;
    }

    return false;
}

template <typename T>
void TriangleSet<T>::readAsciiSTL(FILE *fp, TriangleSelector<T> *selector)
{
    char buf[PLB_CBUFSIZ];
    char *cp;

#ifdef PLB_DEBUG
    char *sp = fgets(buf, PLB_CBUFSIZ, fp);
#else
    (void)fgets(buf, PLB_CBUFSIZ, fp);  // TODO: (void) does not suppress warning. We should alway
                                        // check return. throw exception.
#endif
    PLB_ASSERT(sp != NULL);                   // The input file is badly structured.
    PLB_ASSERT(checkForBufferOverflow(buf));  // Problem with reading one line of text.

    char fmt[32];
    bool failed = false;
    if (sizeof(T) == sizeof(float))
        strcpy(fmt, "%f%f%f");
    else if (sizeof(T) == sizeof(double))
        strcpy(fmt, "%lf%lf%lf");
    else if (sizeof(T) == sizeof(long double))
        strcpy(fmt, "%Lf%Lf%Lf");
    else
        failed = true;

    PLB_ASSERT(!failed);  // The input file cannot be read.

    cp = strstr(buf, "solid");
    PLB_ASSERT(cp != NULL);  // The input file is badly structured.

    plint partId = 0;

    while (cp != NULL && !failed) {
        if (fgets(buf, PLB_CBUFSIZ, fp) == NULL) {
            failed = true;
            PLB_ASSERT(!failed);  // The input file cannot be read.
        }
        PLB_ASSERT(checkForBufferOverflow(buf));  // Problem with reading one line of text.

        if ((cp = strstr(buf, "color")) != NULL) {
            if (fgets(buf, PLB_CBUFSIZ, fp) == NULL) {
                failed = true;
                PLB_ASSERT(!failed);  // The input file cannot be read.
            }
            PLB_ASSERT(checkForBufferOverflow(buf));  // Problem with reading one line of text.
        }

        do {
            if ((cp = strstr(buf, "facet normal")) == NULL) {
                failed = true;
                PLB_ASSERT(!failed);  // The input file cannot be read.
            }
            cp += 12;
            Array<T, 3> n;
            if (sscanf(cp, fmt, &n[0], &n[1], &n[2]) != 3) {
                failed = true;
                PLB_ASSERT(!failed);  // The input file cannot be read.
            }

            if (fgets(buf, PLB_CBUFSIZ, fp) == NULL || strstr(buf, "outer loop") == NULL) {
                failed = true;
                PLB_ASSERT(!failed);  // The input file cannot be read.
            }
            PLB_ASSERT(checkForBufferOverflow(buf));  // Problem with reading one line of text.

            Triangle triangle;
            for (int i = 0; i < 3; i++) {
                if (fgets(buf, PLB_CBUFSIZ, fp) == NULL || (cp = strstr(buf, "vertex")) == NULL) {
                    failed = true;
                    PLB_ASSERT(!failed);  // The input file cannot be read.
                }
                PLB_ASSERT(checkForBufferOverflow(buf));  // Problem with reading one line of text.
                cp += 6;
                triangle[i][0] = T();
                triangle[i][1] = T();
                triangle[i][2] = T();
                if (sscanf(cp, fmt, &triangle[i][0], &triangle[i][1], &triangle[i][2]) != 3) {
                    failed = true;
                    PLB_ASSERT(!failed);  // The input file cannot be read.
                }
            }

            if (fgets(buf, PLB_CBUFSIZ, fp) == NULL || strstr(buf, "endloop") == NULL) {
                failed = true;
                PLB_ASSERT(!failed);  // The input file cannot be read.
            }
            PLB_ASSERT(checkForBufferOverflow(buf));  // Problem with reading one line of text.

            if (fgets(buf, PLB_CBUFSIZ, fp) == NULL || strstr(buf, "endfacet") == NULL) {
                failed = true;
                PLB_ASSERT(!failed);  // The input file cannot be read.
            }
            PLB_ASSERT(checkForBufferOverflow(buf));  // Problem with reading one line of text.

            // We keep triangles that have a zero area, since they might be used to
            // fix topological problems with the connectivity of the triangle set.
            // However, we do not keep triangles that have one or more edges with a
            // zero length. This has as a result, that if a triangle is kept in the
            // triangle set and has zero area, its vertices are distinct but belong
            // to a straight line.
            if (!triangleHasZeroLengthEdges(triangle, eps)) {
                fixOrientation(triangle, n);
                if (selector == 0 || (*selector)(triangle, partId)) {
                    triangles.push_back(triangle);
                }
            }

            if (fgets(buf, PLB_CBUFSIZ, fp) == NULL) {
                failed = true;
                PLB_ASSERT(!failed);  // The input file cannot be read.
            }
            cp = strstr(buf, "endsolid");
            if (cp == NULL) {
                PLB_ASSERT(checkForBufferOverflow(buf));  // Problem with reading one line of text.
            }
        } while (cp == NULL && !failed);

        if (fgets(buf, PLB_CBUFSIZ, fp) == NULL) {
            break;
        }

        cp = strstr(buf, "solid");
        partId++;
    }
}

template <typename T>
void TriangleSet<T>::readBinarySTL(FILE *fp, TriangleSelector<T> *selector)
{
    char buf[PLB_CBUFSIZ];
    std::fill(buf, buf + PLB_CBUFSIZ, '\0');
    unsigned int nt;
    float array[3];
    unsigned short abc;
    bool failed = false;

    int count = 0;
    while (fread(buf, sizeof(char), 80, fp) == 80 && fread(&nt, sizeof(unsigned int), 1, fp) == 1
           && !failed)
    {
        count++;
        plint partId = (plint)count - 1;
        for (unsigned it = 0; it < nt && !failed; it++) {
            if (fread(array, sizeof(float), 3, fp) != 3) {
                failed = true;
                PLB_ASSERT(!failed);  // The input file is badly structured.
            }
            Array<T, 3> n;
            n[0] = array[0];
            n[1] = array[1];
            n[2] = array[2];

            Triangle triangle;
            for (int i = 0; i < 3 && !failed; i++) {
                if (fread(array, sizeof(float), 3, fp) != 3) {
                    failed = true;
                    PLB_ASSERT(!failed);  // The input file is badly structured.
                }
                triangle[i][0] = T();
                triangle[i][1] = T();
                triangle[i][2] = T();
                triangle[i][0] = array[0];
                triangle[i][1] = array[1];
                triangle[i][2] = array[2];
            }

            if (fread(&abc, sizeof(unsigned short), 1, fp) != 1) {
                failed = true;
                PLB_ASSERT(!failed);  // The input file is badly structured.
            }

            // We keep triangles that have a zero area, since they might be used to
            // fix topological problems with the connectivity of the triangle set.
            // However, we do not keep triangles that have one or more edges with a
            // zero length. This has as a result, that if a triangle is kept in the
            // triangle set and has zero area, its vertices are distinct but belong
            // to a straight line.
            if (!triangleHasZeroLengthEdges(triangle, eps)) {
                fixOrientation(triangle, n);
                if (selector == 0 || (*selector)(triangle, partId)) {
                    triangles.push_back(triangle);
                }
            }
        }
    }

    if (count == 0)
        failed = true;

    PLB_ASSERT(!failed);  // The input file is badly structured.
}

template <typename T>
void TriangleSet<T>::readOFF(std::string fname, TriangleSelector<T> *selector)
{
    FILE *fp = fopen(fname.c_str(), "rb");
    PLB_ASSERT(fp != 0);  // The input file cannot be read.

    char buf[PLB_CBUFSIZ];
#ifdef PLB_DEBUG
    char *sp = fgets(buf, PLB_CBUFSIZ, fp);
#else
    (void)fgets(buf, PLB_CBUFSIZ, fp);  // TODO: (void) does not suppress warning. We should alway
                                        // check return. throw exception.
#endif
    PLB_ASSERT(sp != NULL);  // The input file cannot be read.

    char *cp = NULL;

    // Currently only ASCII files with header OFF can be read.
    cp = strstr(buf, "BINARY");
    PLB_ASSERT(cp == NULL);
    std::vector<std::string> prefixes;
    prefixes.push_back("ST");
    prefixes.push_back("C");
    prefixes.push_back("N");
    prefixes.push_back("4");
    prefixes.push_back("n");
    for (int iPrefix = 0; iPrefix < (int)prefixes.size(); iPrefix++) {
        cp = strstr(buf, prefixes[iPrefix].c_str());
        PLB_ASSERT(cp == NULL);
    }

    cp = strstr(buf, "OFF");

    if (cp != NULL) {
        readAsciiOFF(fp, selector);
    } else {
        PLB_ASSERT(false);  // Nothing else is currently supported.
    }

    fclose(fp);
}

template <typename T>
void TriangleSet<T>::readAsciiOFF(FILE *fp, TriangleSelector<T> *selector)
{
    char commentCharacter = '#';
    long NVertices = 0, NFaces = 0, NEdges = 0;
    if (readAhead(fp, commentCharacter) == EOF) {
        PLB_ASSERT(false);  // The input file is badly structured.
    }
    if (fscanf(fp, "%ld%ld%ld", &NVertices, &NFaces, &NEdges) != 3) {
        PLB_ASSERT(false);  // The input file is badly structured.
    }

    char fmt[32];
    if (sizeof(T) == sizeof(float)) {
        strcpy(fmt, "%f%f%f");
    } else if (sizeof(T) == sizeof(double)) {
        strcpy(fmt, "%lf%lf%lf");
    } else if (sizeof(T) == sizeof(long double)) {
        strcpy(fmt, "%Lf%Lf%Lf");
    } else {
        PLB_ASSERT(false);  // The input file cannot be read.
    }

    std::vector<Array<T, 3> > vertices(NVertices);
    for (long iVertex = 0; iVertex < NVertices; iVertex++) {
        if (readAhead(fp, commentCharacter) == EOF) {
            PLB_ASSERT(false);  // The input file is badly structured.
        }
        if (fscanf(fp, fmt, &vertices[iVertex][0], &vertices[iVertex][1], &vertices[iVertex][2])
            != 3) {
            PLB_ASSERT(false);  // The input file is badly structured.
        }
    }

    plint partId = 0;  // Only one part in OFF files.

    triangles.clear();
    for (long iFace = 0; iFace < NFaces; iFace++) {
        if (readAhead(fp, commentCharacter) == EOF) {
            PLB_ASSERT(false);  // The input file is badly structured.
        }
        long nv;
        if (fscanf(fp, "%ld", &nv) != 1) {
            PLB_ASSERT(false);  // The input file is badly structured.
        }
        PLB_ASSERT(nv == 3);  // The surface mesh is not triangulated.
        long ind[3];
        if (fscanf(fp, "%ld%ld%ld", &ind[0], &ind[1], &ind[2]) != 3) {
            PLB_ASSERT(false);  // The input file is badly structured.
        }
        PLB_ASSERT(ind[0] >= 0 && ind[0] < NVertices);
        PLB_ASSERT(ind[1] >= 0 && ind[1] < NVertices);
        PLB_ASSERT(ind[2] >= 0 && ind[2] < NVertices);

        Triangle triangle;
        triangle[0] = vertices[ind[0]];
        triangle[1] = vertices[ind[1]];
        triangle[2] = vertices[ind[2]];

        // We keep triangles that have a zero area, since they might be used to
        // fix topological problems with the connectivity of the triangle set.
        // However, we do not keep triangles that have one or more edges with a
        // zero length. This has as a result, that if a triangle is kept in the
        // triangle set and has zero area, its vertices are distinct but belong
        // to a straight line.
        if (!triangleHasZeroLengthEdges(triangle, eps)) {
            if (selector == 0 || (*selector)(triangle, partId)) {
                triangles.push_back(triangle);
            }
        }
    }
}

template <typename T>
bool TriangleSet<T>::triangleHasZeroArea(Triangle const &triangle, T epsilon) const
{
    T area = computeTriangleArea(triangle[0], triangle[1], triangle[2]);
    if (util::isZero(area, epsilon)) {
        return true;
    }
    return false;
}

template <typename T>
bool TriangleSet<T>::triangleHasZeroLengthEdges(Triangle const &triangle, T epsilon) const
{
    if (util::isZero(norm(triangle[1] - triangle[0]), epsilon)
        || util::isZero(norm(triangle[2] - triangle[0]), epsilon)
        || util::isZero(norm(triangle[2] - triangle[1]), epsilon))
    {
        return true;
    }
    return false;
}

template <typename T>
void TriangleSet<T>::fixOrientation(Triangle &triangle, Array<T, 3> const &n) const
{
    bool isAreaWeighted = false;
    Array<T, 3> computedNormal =
        computeTriangleNormal(triangle[0], triangle[1], triangle[2], isAreaWeighted);
    if (dot(computedNormal, n) < (T)0) {
        std::swap(triangle[1], triangle[2]);
    }
}

template <typename T>
bool TriangleSet<T>::areTheSameTriangle(Triangle const &t1, Triangle const &t2, T epsilon) const
{
    Array<T, 3> const &v10 = t1[0];
    Array<T, 3> const &v11 = t1[1];
    Array<T, 3> const &v12 = t1[2];

    Array<T, 3> const &v20 = t2[0];
    Array<T, 3> const &v21 = t2[1];
    Array<T, 3> const &v22 = t2[2];

    T n1020 = norm(v10 - v20);
    T n1021 = norm(v10 - v21);
    T n1022 = norm(v10 - v22);
    if (!util::isZero(n1020, epsilon) && !util::isZero(n1021, epsilon)
        && !util::isZero(n1022, epsilon)) {
        return false;
    }

    T n1120 = norm(v11 - v20);
    T n1121 = norm(v11 - v21);
    T n1122 = norm(v11 - v22);
    if (!util::isZero(n1120, epsilon) && !util::isZero(n1121, epsilon)
        && !util::isZero(n1122, epsilon)) {
        return false;
    }

    T n1220 = norm(v12 - v20);
    T n1221 = norm(v12 - v21);
    T n1222 = norm(v12 - v22);
    if (!util::isZero(n1220, epsilon) && !util::isZero(n1221, epsilon)
        && !util::isZero(n1222, epsilon)) {
        return false;
    }

    return true;
}

template <typename T>
bool TriangleSet<T>::containedAndErase(
    Triangle const &triangle, std::vector<Triangle> &triangles, T epsilon) const
{
    typename std::vector<Triangle>::iterator it;
    for (it = triangles.begin(); it != triangles.end(); ++it) {
        if (areTheSameTriangle(triangle, *it, epsilon)) {
            triangles.erase(it);
            return true;
        }
    }
    return false;
}

template <typename T>
void TriangleSet<T>::translate(Array<T, 3> const &vector)
{
    if (util::isZero(norm(vector)))
        return;

    plint size = triangles.size();
    if (size == 0)
        return;

    for (plint i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            triangles[i][j] += vector;
        }
    }

    boundingCuboid.lowerLeftCorner += vector;
    boundingCuboid.upperRightCorner += vector;
}

template <typename T>
void TriangleSet<T>::scale(T alpha)
{
    if (util::isOne(alpha))
        return;

    plint size = triangles.size();
    if (size == 0)
        return;

    for (plint i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            triangles[i][j] *= alpha;
        }
    }

    computeMinMaxEdges();
    computeMinMaxAreas();
    computeBoundingCuboid();
}

template <typename T>
void TriangleSet<T>::scale(T alpha, T beta, T gamma)
{
    if (util::isOne(alpha) && util::isOne(beta) && util::isOne(gamma))
        return;

    plint size = triangles.size();
    if (size == 0)
        return;

    for (plint i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            triangles[i][j][0] *= alpha;
            triangles[i][j][1] *= beta;
            triangles[i][j][2] *= gamma;
        }
    }

    computeMinMaxEdges();
    computeMinMaxAreas();
    computeBoundingCuboid();
}

template <typename T>
void TriangleSet<T>::rotateAtOrigin(Array<T, 3> const &normedAxis, T theta)
{
    plint size = (plint)triangles.size();
    for (plint iTriangle = 0; iTriangle < size; iTriangle++) {
        for (int iVertex = 0; iVertex < 3; iVertex++) {
            triangles[iTriangle][iVertex] =
                plb::rotateAtOrigin(triangles[iTriangle][iVertex], normedAxis, theta);
        }
    }

    computeBoundingCuboid();
}

template <typename T>
void TriangleSet<T>::rotateAxesTo(Array<T, 3> const &ex, Array<T, 3> const &ey)
{
    Array<T, 3> ez(crossProduct(ex, ey));
    plint size = (plint)triangles.size();
    for (plint iTriangle = 0; iTriangle < size; iTriangle++) {
        for (int iVertex = 0; iVertex < 3; iVertex++) {
            Array<T, 3> &v = triangles[iTriangle][iVertex];
            v = Array<T, 3>(
                ex[0] * v[0] + ey[0] * v[1] + ez[0] * v[2],
                ex[1] * v[0] + ey[1] * v[1] + ez[1] * v[2],
                ex[2] * v[0] + ey[2] * v[1] + ez[2] * v[2]);
        }
    }

    computeBoundingCuboid();
}

template <typename T>
void TriangleSet<T>::rotateAxesFrom(Array<T, 3> const &ex, Array<T, 3> const &ey)
{
    Array<T, 3> ez(crossProduct(ex, ey));
    plint size = (plint)triangles.size();
    for (plint iTriangle = 0; iTriangle < size; iTriangle++) {
        for (int iVertex = 0; iVertex < 3; iVertex++) {
            Array<T, 3> &v = triangles[iTriangle][iVertex];
            v = Array<T, 3>(
                ex[0] * v[0] + ex[1] * v[1] + ex[2] * v[2],
                ey[0] * v[0] + ey[1] * v[1] + ey[2] * v[2],
                ez[0] * v[0] + ez[1] * v[1] + ez[2] * v[2]);
        }
    }

    computeBoundingCuboid();
}

template <typename T>
void TriangleSet<T>::rotate(T phi, T theta, T psi)
{
#ifdef PLB_DEBUG
    static T pi = std::acos((T)-1);
#endif
    PLB_ASSERT(util::greaterEqual(theta, (T)0, eps) && util::lessEqual(theta, pi, eps));

    plint size = triangles.size();
    if (size == 0)
        return;

    T cosPhi = std::cos(phi);
    T sinPhi = std::sin(phi);
    T cosTheta = std::cos(theta);
    T sinTheta = std::sin(theta);
    T cosPsi = std::cos(psi);
    T sinPsi = std::sin(psi);

    T a[3][3];
    a[0][0] = (T)1.0;
    a[0][1] = (T)0.0;
    a[0][2] = (T)0.0;
    a[1][0] = (T)0.0;
    a[1][1] = cosTheta;
    a[1][2] = -sinTheta;
    a[2][0] = (T)0.0;
    a[2][1] = sinTheta;
    a[2][2] = cosTheta;

    T b[3][3];
    b[0][0] = cosPhi;
    b[0][1] = -sinPhi;
    b[0][2] = (T)0.0;
    b[1][0] = sinPhi;
    b[1][1] = cosPhi;
    b[1][2] = (T)0.0;
    b[2][0] = (T)0.0;
    b[2][1] = (T)0.0;
    b[2][2] = (T)1.0;

    T c[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            c[i][j] = (T)0.0;
            for (int k = 0; k < 3; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    b[0][0] = cosPsi;
    b[0][1] = -sinPsi;
    b[0][2] = (T)0.0;
    b[1][0] = sinPsi;
    b[1][1] = cosPsi;
    b[1][2] = (T)0.0;
    b[2][0] = (T)0.0;
    b[2][1] = (T)0.0;
    b[2][2] = (T)1.0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i][j] = (T)0.0;
            for (int k = 0; k < 3; k++) {
                a[i][j] += b[i][k] * c[k][j];
            }
        }
    }

    for (plint iTriangle = 0; iTriangle < size; iTriangle++) {
        for (int iVertex = 0; iVertex < 3; iVertex++) {
            Array<T, 3> x = triangles[iTriangle][iVertex];
            for (int i = 0; i < 3; i++) {
                triangles[iTriangle][iVertex][i] = (T)0.0;
                for (int j = 0; j < 3; j++) {
                    triangles[iTriangle][iVertex][i] += a[i][j] * x[j];
                }
            }
        }
    }

    computeBoundingCuboid();
}

template <typename T>
void TriangleSet<T>::clear()
{
    std::vector<Triangle>().swap(triangles);
    minEdgeLength = std::numeric_limits<T>::max();
    maxEdgeLength = std::numeric_limits<T>::min();
    minTriangleArea = std::numeric_limits<T>::max();
    maxTriangleArea = std::numeric_limits<T>::min();

    T huge = std::numeric_limits<T>::max();

    boundingCuboid.lowerLeftCorner = Array<T, 3>((T)huge, (T)huge, (T)huge);
    boundingCuboid.upperRightCorner = Array<T, 3>((T)-huge, (T)-huge, (T)-huge);

    eps = (T)0.0;
}

template <typename T>
void TriangleSet<T>::merge(std::vector<TriangleSet<T> *> meshes)
{
    PLB_ASSERT(meshes.size() != 0);

    triangles.assign(meshes[0]->getTriangles().begin(), meshes[0]->getTriangles().end());
    minEdgeLength = meshes[0]->getMinEdgeLength();
    maxEdgeLength = meshes[0]->getMaxEdgeLength();
    minTriangleArea = meshes[0]->getMinTriangleArea();
    maxTriangleArea = meshes[0]->getMaxTriangleArea();
    boundingCuboid = meshes[0]->getBoundingCuboid();
    for (pluint i = 1; i < meshes.size(); i++) {
        triangles.insert(
            triangles.end(), meshes[i]->getTriangles().begin(), meshes[i]->getTriangles().end());
        minEdgeLength = std::min(minEdgeLength, meshes[i]->getMinEdgeLength());
        maxEdgeLength = std::max(maxEdgeLength, meshes[i]->getMaxEdgeLength());
        minTriangleArea = std::min(minTriangleArea, meshes[i]->getMinTriangleArea());
        maxTriangleArea = std::max(maxTriangleArea, meshes[i]->getMaxTriangleArea());

        Cuboid<T> bcuboid = meshes[i]->getBoundingCuboid();
        for (plint j = 0; j < 3; j++) {
            boundingCuboid.lowerLeftCorner[j] =
                std::min(boundingCuboid.lowerLeftCorner[j], bcuboid.lowerLeftCorner[j]);
            boundingCuboid.upperRightCorner[j] =
                std::max(boundingCuboid.upperRightCorner[j], bcuboid.upperRightCorner[j]);
        }
    }
}

template <typename T>
void TriangleSet<T>::append(TriangleSet<T> const &mesh)
{
    triangles.insert(triangles.end(), mesh.getTriangles().begin(), mesh.getTriangles().end());
    minEdgeLength = std::min(minEdgeLength, mesh.getMinEdgeLength());
    maxEdgeLength = std::max(maxEdgeLength, mesh.getMaxEdgeLength());
    minTriangleArea = std::min(minTriangleArea, mesh.getMinTriangleArea());
    maxTriangleArea = std::max(maxTriangleArea, mesh.getMaxTriangleArea());

    Cuboid<T> bcuboid = mesh.getBoundingCuboid();
    for (plint j = 0; j < 3; j++) {
        boundingCuboid.lowerLeftCorner[j] =
            std::min(boundingCuboid.lowerLeftCorner[j], bcuboid.lowerLeftCorner[j]);
        boundingCuboid.upperRightCorner[j] =
            std::max(boundingCuboid.upperRightCorner[j], bcuboid.upperRightCorner[j]);
    }
}

template <typename T>
void TriangleSet<T>::refine()
{
    std::vector<Triangle> newTriangles;
    std::vector<Triangle> newZeroAreaTriangles;

    for (pluint i = 0; i < triangles.size(); ++i) {
        Triangle const &triangle = triangles[i];

        Array<T, 3> v00 = triangle[0];
        Array<T, 3> v11 = triangle[1];
        Array<T, 3> v22 = triangle[2];

        Array<T, 3> v01 = (T)0.5 * (v00 + v11);
        Array<T, 3> v12 = (T)0.5 * (v11 + v22);
        Array<T, 3> v20 = (T)0.5 * (v22 + v00);

        Triangle newTriangle;

        newTriangle[0] = v01;
        newTriangle[1] = v12;
        newTriangle[2] = v20;
        if (triangleHasZeroArea(triangle, eps)) {
            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
            {
                newZeroAreaTriangles.push_back(newTriangle);
            }
        } else {
            newTriangles.push_back(newTriangle);
        }

        newTriangle[0] = v00;
        newTriangle[1] = v01;
        newTriangle[2] = v20;
        if (triangleHasZeroArea(triangle, eps)) {
            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
            {
                newZeroAreaTriangles.push_back(newTriangle);
            }
        } else {
            newTriangles.push_back(newTriangle);
        }

        newTriangle[0] = v01;
        newTriangle[1] = v11;
        newTriangle[2] = v12;
        if (triangleHasZeroArea(triangle, eps)) {
            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
            {
                newZeroAreaTriangles.push_back(newTriangle);
            }
        } else {
            newTriangles.push_back(newTriangle);
        }

        newTriangle[0] = v20;
        newTriangle[1] = v12;
        newTriangle[2] = v22;
        if (triangleHasZeroArea(triangle, eps)) {
            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
            {
                newZeroAreaTriangles.push_back(newTriangle);
            }
        } else {
            newTriangles.push_back(newTriangle);
        }
    }

    newTriangles.insert(
        newTriangles.end(), newZeroAreaTriangles.begin(), newZeroAreaTriangles.end());

    triangles.swap(newTriangles);

    computeMinMaxEdges();
    computeMinMaxAreas();
    // Don't recompute bounding cuboid, because it is unchanged.
}

template <typename T>
void TriangleSet<T>::refine(T edgeLengthThreshold)
{
    PLB_ASSERT(util::greaterThan(edgeLengthThreshold, (T)0, eps));

    T thr2 = util::sqr(edgeLengthThreshold);

    std::vector<Triangle> newTriangles;
    std::vector<Triangle> newZeroAreaTriangles;
    Triangle newTriangle;

    for (pluint i = 0; i < triangles.size(); ++i) {
        Triangle const &triangle = triangles[i];

        Array<T, 3> v0 = triangle[0];
        Array<T, 3> v1 = triangle[1];
        Array<T, 3> v2 = triangle[2];

        T e0 = normSqr<T, 3>(v1 - v0);
        T e1 = normSqr<T, 3>(v2 - v1);
        T e2 = normSqr<T, 3>(v0 - v2);

        if (e0 >= thr2) {
            Array<T, 3> vA = (T)0.5 * (v0 + v1);
            if (e1 >= thr2) {
                Array<T, 3> vB = (T)0.5 * (v1 + v2);
                if (e2 >= thr2) {
                    Array<T, 3> vC = (T)0.5 * (v2 + v0);

                    newTriangle[0] = v0;
                    newTriangle[1] = vA;
                    newTriangle[2] = vC;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }

                    newTriangle[0] = vA;
                    newTriangle[1] = v1;
                    newTriangle[2] = vB;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }

                    newTriangle[0] = vC;
                    newTriangle[1] = vB;
                    newTriangle[2] = v2;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }

                    newTriangle[0] = vA;
                    newTriangle[1] = vB;
                    newTriangle[2] = vC;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }
                } else {
                    newTriangle[0] = vA;
                    newTriangle[1] = v1;
                    newTriangle[2] = vB;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }

                    T eA2 = normSqr<T, 3>(vA - v2);
                    T eB0 = normSqr<T, 3>(vB - v0);
                    if (eA2 < eB0) {
                        newTriangle[0] = v2;
                        newTriangle[1] = v0;
                        newTriangle[2] = vA;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = v2;
                        newTriangle[1] = vA;
                        newTriangle[2] = vB;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangle[0] = v0;
                        newTriangle[1] = vA;
                        newTriangle[2] = vB;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = v0;
                        newTriangle[1] = vB;
                        newTriangle[2] = v2;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    }
                }
            } else {
                if (e2 >= thr2) {
                    Array<T, 3> vC = (T)0.5 * (v2 + v0);

                    newTriangle[0] = v0;
                    newTriangle[1] = vA;
                    newTriangle[2] = vC;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }

                    T eA2 = normSqr<T, 3>(vA - v2);
                    T eC1 = normSqr<T, 3>(vC - v1);
                    if (eA2 < eC1) {
                        newTriangle[0] = v2;
                        newTriangle[1] = vC;
                        newTriangle[2] = vA;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = v2;
                        newTriangle[1] = vA;
                        newTriangle[2] = v1;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangle[0] = v1;
                        newTriangle[1] = v2;
                        newTriangle[2] = vC;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = v1;
                        newTriangle[1] = vC;
                        newTriangle[2] = vA;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    }
                } else {
                    newTriangle[0] = v2;
                    newTriangle[1] = v0;
                    newTriangle[2] = vA;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }

                    newTriangle[0] = v2;
                    newTriangle[1] = vA;
                    newTriangle[2] = v1;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }
                }
            }
        } else {
            if (e1 >= thr2) {
                Array<T, 3> vB = (T)0.5 * (v1 + v2);
                if (e2 >= thr2) {
                    Array<T, 3> vC = (T)0.5 * (v2 + v0);

                    newTriangle[0] = v2;
                    newTriangle[1] = vC;
                    newTriangle[2] = vB;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }

                    T eB0 = normSqr<T, 3>(vB - v0);
                    T eC1 = normSqr<T, 3>(vC - v1);
                    if (eB0 < eC1) {
                        newTriangle[0] = v0;
                        newTriangle[1] = v1;
                        newTriangle[2] = vB;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = v0;
                        newTriangle[1] = vB;
                        newTriangle[2] = vC;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangle[0] = v1;
                        newTriangle[1] = vC;
                        newTriangle[2] = v0;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }

                        newTriangle[0] = v1;
                        newTriangle[1] = vB;
                        newTriangle[2] = vC;
                        if (triangleHasZeroArea(triangle, eps)) {
                            if (!triangleHasZeroLengthEdges(newTriangle, eps)
                                && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                            {
                                newZeroAreaTriangles.push_back(newTriangle);
                            }
                        } else {
                            newTriangles.push_back(newTriangle);
                        }
                    }
                } else {
                    newTriangle[0] = v0;
                    newTriangle[1] = v1;
                    newTriangle[2] = vB;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }

                    newTriangle[0] = v0;
                    newTriangle[1] = vB;
                    newTriangle[2] = v2;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }
                }
            } else {
                if (e2 >= thr2) {
                    Array<T, 3> vC = (T)0.5 * (v2 + v0);

                    newTriangle[0] = v1;
                    newTriangle[1] = vC;
                    newTriangle[2] = v0;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }

                    newTriangle[0] = v1;
                    newTriangle[1] = v2;
                    newTriangle[2] = vC;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }
                } else {
                    newTriangle[0] = v0;
                    newTriangle[1] = v1;
                    newTriangle[2] = v2;
                    if (triangleHasZeroArea(triangle, eps)) {
                        if (!triangleHasZeroLengthEdges(newTriangle, eps)
                            && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                        {
                            newZeroAreaTriangles.push_back(newTriangle);
                        }
                    } else {
                        newTriangles.push_back(newTriangle);
                    }
                }
            }
        }
    }

    newTriangles.insert(
        newTriangles.end(), newZeroAreaTriangles.begin(), newZeroAreaTriangles.end());

    triangles.swap(newTriangles);

    computeMinMaxEdges();
    computeMinMaxAreas();
    // Don't recompute bounding cuboid, because it is unchanged.
}

template <typename T>
bool TriangleSet<T>::refineRecursively(T targetMaxEdgeLength, plint maxNumIterations)
{
    PLB_ASSERT(util::greaterThan(targetMaxEdgeLength, (T)0, eps));
    PLB_ASSERT(maxNumIterations > 0);

    plint iter = 0;
    while (maxEdgeLength >= targetMaxEdgeLength && iter < maxNumIterations) {
        refine(targetMaxEdgeLength);
        iter++;
    }

    if (maxEdgeLength < targetMaxEdgeLength) {
        return true;
    }

    return false;
}

template <typename T>
void TriangleSet<T>::refineByArea(T triangleAreaThreshold)
{
    PLB_ASSERT(util::greaterThan(triangleAreaThreshold, (T)0, eps));

    std::vector<Triangle> newTriangles;
    std::vector<Triangle> newZeroAreaTriangles;

    for (pluint i = 0; i < triangles.size(); ++i) {
        Triangle const &triangle = triangles[i];
        T area = computeTriangleArea(triangle[0], triangle[1], triangle[2]);
        if (area >= triangleAreaThreshold) {
            Array<T, 3> v00 = triangle[0];
            Array<T, 3> v11 = triangle[1];
            Array<T, 3> v22 = triangle[2];

            Array<T, 3> v01 = (T)0.5 * (v00 + v11);
            Array<T, 3> v12 = (T)0.5 * (v11 + v22);
            Array<T, 3> v20 = (T)0.5 * (v22 + v00);

            Triangle newTriangle;

            newTriangle[0] = v01;
            newTriangle[1] = v12;
            newTriangle[2] = v20;
            if (triangleHasZeroArea(triangle, eps)) {
                if (!triangleHasZeroLengthEdges(newTriangle, eps)
                    && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                {
                    newZeroAreaTriangles.push_back(newTriangle);
                }
            } else {
                newTriangles.push_back(newTriangle);
            }

            newTriangle[0] = v00;
            newTriangle[1] = v01;
            newTriangle[2] = v20;
            if (triangleHasZeroArea(triangle, eps)) {
                if (!triangleHasZeroLengthEdges(newTriangle, eps)
                    && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                {
                    newZeroAreaTriangles.push_back(newTriangle);
                }
            } else {
                newTriangles.push_back(newTriangle);
            }

            newTriangle[0] = v01;
            newTriangle[1] = v11;
            newTriangle[2] = v12;
            if (triangleHasZeroArea(triangle, eps)) {
                if (!triangleHasZeroLengthEdges(newTriangle, eps)
                    && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                {
                    newZeroAreaTriangles.push_back(newTriangle);
                }
            } else {
                newTriangles.push_back(newTriangle);
            }

            newTriangle[0] = v20;
            newTriangle[1] = v12;
            newTriangle[2] = v22;
            if (triangleHasZeroArea(triangle, eps)) {
                if (!triangleHasZeroLengthEdges(newTriangle, eps)
                    && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                {
                    newZeroAreaTriangles.push_back(newTriangle);
                }
            } else {
                newTriangles.push_back(newTriangle);
            }
        } else {
            newTriangles.push_back(triangle);
        }
    }

    newTriangles.insert(
        newTriangles.end(), newZeroAreaTriangles.begin(), newZeroAreaTriangles.end());

    triangles.swap(newTriangles);

    computeMinMaxEdges();
    computeMinMaxAreas();
    // Don't recompute bounding cuboid, because it is unchanged.
}

template <typename T>
bool TriangleSet<T>::refineByAreaRecursively(T targetMaxTriangleArea, plint maxNumIterations)
{
    PLB_ASSERT(util::greaterThan(targetMaxTriangleArea, (T)0, eps));
    PLB_ASSERT(maxNumIterations > 0);

    plint iter = 0;
    while (maxTriangleArea >= targetMaxTriangleArea && iter < maxNumIterations) {
        refineByArea(targetMaxTriangleArea);
        iter++;
    }

    if (maxTriangleArea < targetMaxTriangleArea) {
        return true;
    }

    return false;
}

template <typename T>
void TriangleSet<T>::select(TriangleSelector<T> const &selector)
{
    plint partId = 0;  // There is only one part in the TriangleSet.
    std::vector<Triangle> newTriangles;

    for (pluint i = 0; i < triangles.size(); ++i) {
        Triangle const &triangle = triangles[i];
        if (selector(triangle, partId)) {
            newTriangles.push_back(triangle);
        }
    }

    triangles.swap(newTriangles);

    computeMinMaxEdges();
    computeMinMaxAreas();
    computeBoundingCuboid();
}

template <typename T>
void TriangleSet<T>::removeTrianglesWithOrientation(Array<T, 3> const &normal, T tolerance)
{
    std::vector<Triangle> newTriangles;
    for (pluint i = 0; i < triangles.size(); ++i) {
        Triangle triangle = triangles[i];
        Array<T, 3> triangleNormal =
            computeTriangleNormal(triangle[0], triangle[1], triangle[2], false);
        if (dot(triangleNormal, normal) < tolerance) {  // We decide to keep zero-area triangles.
            newTriangles.push_back(triangle);
        }
    }
    triangles.swap(newTriangles);
    computeMinMaxEdges();
    computeMinMaxAreas();
    computeBoundingCuboid();
}

template <typename T>
void TriangleSet<T>::reverseOrientation()
{
    plint size = triangles.size();
    for (plint i = 0; i < size; i++)
        std::swap(triangles[i][1], triangles[i][2]);
}

template <typename T>
void TriangleSet<T>::toFloatAndBack()
{
    plint size = triangles.size();
    for (plint i = 0; i < size; i++) {
        Triangle &triangle = triangles[i];
        for (plint iVertex = 0; iVertex < 3; iVertex++) {
            Array<T, 3> &vertex = triangle[iVertex];
            vertex[0] = (T)(float)vertex[0];
            vertex[1] = (T)(float)vertex[1];
            vertex[2] = (T)(float)vertex[2];
        }
    }
}

template <typename T>
bool TriangleSet<T>::hasZeroAreaTriangles() const
{
    plint size = triangles.size();
    for (plint i = 0; i < size; i++) {
        Triangle const &triangle = triangles[i];
        if (triangleHasZeroArea(triangle, eps)) {
            return true;
        }
    }
    return false;
}

template <typename T>
plint TriangleSet<T>::numZeroAreaTriangles() const
{
    plint size = triangles.size();
    plint count = 0;
    for (plint i = 0; i < size; i++) {
        Triangle const &triangle = triangles[i];
        if (triangleHasZeroArea(triangle, eps)) {
            count++;
        }
    }
    return count;
}

template <typename T>
bool TriangleSet<T>::hasFloatingPointPrecisionDependence() const
{
    if (sizeof(T) == sizeof(float)) {
        return false;  // The precision used is the lowest, so the check cannot be made.
    }

    // Here we use epsFLT = eps, because eps (the precision) is chosen by the user
    // independently of the floating point precision of the triangle set (from the
    // STL file). As an example, the user may have chosen a triangle set from a
    // binary STL file (which contains floats), have the code compiled with long doubles,
    // and provide DBL as a precision to determine eps. All checks should be made
    // with the precision (eps) explicitly chosen by the user.
    T epsFLT = eps;
    plint size = triangles.size();
    for (plint i = 0; i < size; i++) {
        Triangle const &triangle = triangles[i];

        Triangle triangleFLT = triangle;
        for (plint iVertex = 0; iVertex < 3; iVertex++) {
            Array<T, 3> &vertex = triangleFLT[iVertex];
            vertex[0] = (T)(float)vertex[0];
            vertex[1] = (T)(float)vertex[1];
            vertex[2] = (T)(float)vertex[2];
        }

        bool rv = triangleHasZeroArea(triangle, eps);
        bool rvFLT = triangleHasZeroArea(triangleFLT, epsFLT);

        if (!util::boolIsEqual(rv, rvFLT)) {
            return true;
        }
    }

    return false;
}

template <typename T>
void TriangleSet<T>::writeAsciiSTL(
    std::string fname, int numDecimalDigits, bool mainProcOnly, TriangleSelector<T> *selector) const
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        if (triangles.empty()) {
            delete selector;
            return;
        }
        FILE *fp = fopen(fname.c_str(), "w");
        PLB_ASSERT(fp != 0);
        T scale = (T)1;
        Array<T, 3> offset(Array<T, 3>::zero());
        writeAsciiContentSTL(fp, "plb", scale, offset, numDecimalDigits, selector);
        fclose(fp);
    } else {
        delete selector;
    }
}

template <typename T>
void TriangleSet<T>::writeAsciiSTL(
    std::string fname, T scale, Array<T, 3> const &offset, int numDecimalDigits, bool mainProcOnly,
    TriangleSelector<T> *selector) const
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        if (triangles.empty()) {
            delete selector;
            return;
        }
        FILE *fp = fopen(fname.c_str(), "w");
        PLB_ASSERT(fp != 0);
        writeAsciiContentSTL(fp, "plb", scale, offset, numDecimalDigits, selector);
        fclose(fp);
    } else {
        delete selector;
    }
}

template <typename T>
void TriangleSet<T>::writeAsciiContentSTL(
    FILE *fp, std::string solidName, T scale, Array<T, 3> const &offset, int numDecimalDigits,
    TriangleSelector<T> *selector) const
{
    plint partId = 0;
    if (selector != 0) {
        plint numSelectedTriangles = 0;
        for (pluint i = 0; i < triangles.size(); i++) {
            if ((*selector)(triangles[i], partId)) {
                numSelectedTriangles++;
            }
        }
        if (numSelectedTriangles == 0) {
            delete selector;
            return;
        }
    }

    std::string header(std::string("solid ") + solidName + std::string("\n"));
    std::string footer(std::string("endsolid ") + solidName + std::string("\n"));

    fprintf(fp, "%s", header.c_str());
    bool isAreaWeighted = false;
    int d = numDecimalDigits;
    for (pluint i = 0; i < triangles.size(); i++) {
        if (selector == 0 || (*selector)(triangles[i], partId)) {
            Array<double, 3> v0 = scale * triangles[i][0] + offset;
            Array<double, 3> v1 = scale * triangles[i][1] + offset;
            Array<double, 3> v2 = scale * triangles[i][2] + offset;

            Array<double, 3> n = computeTriangleNormal(v0, v1, v2, isAreaWeighted);

            fprintf(fp, "  facet normal % .*e % .*e % .*e\n", d, n[0], d, n[1], d, n[2]);
            fprintf(fp, "    outer loop\n");
            fprintf(fp, "      vertex % .*e % .*e % .*e\n", d, v0[0], d, v0[1], d, v0[2]);
            fprintf(fp, "      vertex % .*e % .*e % .*e\n", d, v1[0], d, v1[1], d, v1[2]);
            fprintf(fp, "      vertex % .*e % .*e % .*e\n", d, v2[0], d, v2[1], d, v2[2]);
            fprintf(fp, "    endloop\n");
            fprintf(fp, "  endfacet\n");
        }
    }
    fprintf(fp, "%s", footer.c_str());

    delete selector;
}

template <typename T>
void TriangleSet<T>::writeBinarySTL(
    std::string fname, bool mainProcOnly, TriangleSelector<T> *selector) const
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        if (triangles.empty()) {
            delete selector;
            return;
        }
        FILE *fp = fopen(fname.c_str(), "wb");
        PLB_ASSERT(fp != 0);
        T scale = (T)1;
        Array<T, 3> offset(Array<T, 3>::zero());
        writeBinaryContentSTL(fp, "plb", scale, offset, selector);
        fclose(fp);
    } else {
        delete selector;
    }
}

template <typename T>
void TriangleSet<T>::writeBinarySTL(
    std::string fname, T scale, Array<T, 3> const &offset, bool mainProcOnly,
    TriangleSelector<T> *selector) const
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        if (triangles.empty()) {
            delete selector;
            return;
        }
        FILE *fp = fopen(fname.c_str(), "wb");
        PLB_ASSERT(fp != 0);
        writeBinaryContentSTL(fp, "plb", scale, offset, selector);
        fclose(fp);
    } else {
        delete selector;
    }
}

template <typename T>
void TriangleSet<T>::writeBinaryContentSTL(
    FILE *fp, std::string solidName, T scale, Array<T, 3> const &offset,
    TriangleSelector<T> *selector) const
{
    plint partId = 0;
    if (selector != 0) {
        plint numSelectedTriangles = 0;
        for (pluint i = 0; i < triangles.size(); i++) {
            if ((*selector)(triangles[i], partId)) {
                numSelectedTriangles++;
            }
        }
        if (numSelectedTriangles == 0) {
            delete selector;
            return;
        }
    }

    unsigned int nt = triangles.size();
    unsigned short abc = 0;
    char buf[80] = {'\0'};

    char *name = (char *)malloc((solidName.length() + 1) * sizeof(char));
    strcpy(name, solidName.c_str());
    if (strlen(name) > 79) {
        name[79] = '\0';
    }
    strcpy(buf, name);
    free(name);

    fwrite(buf, sizeof(char), 80, fp);
    fwrite(&nt, sizeof(unsigned int), 1, fp);
    bool isAreaWeighted = false;
    for (unsigned int i = 0; i < nt; i++) {
        if (selector == 0 || (*selector)(triangles[i], partId)) {
            Array<T, 3> v0 = scale * triangles[i][0] + offset;
            Array<T, 3> v1 = scale * triangles[i][1] + offset;
            Array<T, 3> v2 = scale * triangles[i][2] + offset;

            Array<T, 3> nrml = computeTriangleNormal(v0, v1, v2, isAreaWeighted);

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

    delete selector;
}

template <typename T>
int TriangleSet<T>::cutTriangleWithPlane(
    Plane<T> const &plane, Triangle const &triangle, std::vector<Triangle> &newTriangles,
    std::vector<Triangle> &newZeroAreaTriangles) const
{
    int vertexTags[3];

    // Tag the triangle vertices.
    for (int iVertex = 0; iVertex < 3; iVertex++) {
        Array<T, 3> tmp = triangle[iVertex] - plane.point;
        T norm_tmp = norm(tmp);
        if (!util::isZero(norm_tmp, eps)) {
            tmp /= norm_tmp;
        } else {
            tmp[0] = tmp[1] = tmp[2] = (T)0.0;
        }
        T dotp = dot(tmp, plane.normal);
        if (util::isZero(dotp, eps)) {
            vertexTags[iVertex] = 0;
        } else if (util::greaterThan(dotp, (T)0, eps)) {
            vertexTags[iVertex] = -1;
        } else if (util::lessThan(dotp, (T)0, eps)) {
            vertexTags[iVertex] = 1;
        } else {
            return -1;
        }
    }

    // All three vertices belong to one side of the cut plane.
    if (vertexTags[0] == 1 && vertexTags[1] == 1 && vertexTags[2] == 1) {
        if (triangleHasZeroArea(triangle, eps)) {
            if (!triangleHasZeroLengthEdges(triangle, eps)
                && !containedAndErase(triangle, newZeroAreaTriangles, eps))
            {
                newZeroAreaTriangles.push_back(triangle);
            }
        } else {
            newTriangles.push_back(triangle);
        }
        return 1;
    } else if (vertexTags[0] == -1 && vertexTags[1] == -1 && vertexTags[2] == -1) {
        return 0;
    }

    // One vertex belongs to one side of the cut plane and the other two vertices
    //   belong to the other side.
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) != 3 ? (i + 1) : 0;
        int k = (j + 1) != 3 ? (j + 1) : 0;

        if (vertexTags[i] == 1 && vertexTags[j] == -1 && vertexTags[k] == -1) {
            Array<T, 3> intersection_ij((T)0.0, (T)0.0, (T)0.0),
                intersection_ik((T)0.0, (T)0.0, (T)0.0);
            int rv = 0;
            rv =
                lineIntersectionWithPlane<T>(plane, triangle[i], triangle[j], eps, intersection_ij);
            if (rv != 1) {
                return -1;
            }
            rv =
                lineIntersectionWithPlane<T>(plane, triangle[i], triangle[k], eps, intersection_ik);
            if (rv != 1) {
                return -1;
            }
            Triangle newTriangle(triangle[i], intersection_ij, intersection_ik);
            if (triangleHasZeroArea(triangle, eps)) {
                if (!triangleHasZeroLengthEdges(newTriangle, eps)
                    && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                {
                    newZeroAreaTriangles.push_back(newTriangle);
                }
            } else {
                newTriangles.push_back(newTriangle);
            }
            return 1;
        } else if (vertexTags[i] == -1 && vertexTags[j] == 1 && vertexTags[k] == 1) {
            Array<T, 3> intersection_ij((T)0.0, (T)0.0, (T)0.0),
                intersection_ik((T)0.0, (T)0.0, (T)0.0);
            int rv = 0;
            rv =
                lineIntersectionWithPlane<T>(plane, triangle[i], triangle[j], eps, intersection_ij);
            if (rv != 1) {
                return -1;
            }
            rv =
                lineIntersectionWithPlane<T>(plane, triangle[i], triangle[k], eps, intersection_ik);
            if (rv != 1) {
                return -1;
            }
            Triangle newTriangle_0(triangle[k], intersection_ij, triangle[j]);
            Triangle newTriangle_1(triangle[k], intersection_ik, intersection_ij);
            if (triangleHasZeroArea(triangle, eps)) {
                if (!triangleHasZeroLengthEdges(newTriangle_0, eps)
                    && !containedAndErase(newTriangle_0, newZeroAreaTriangles, eps))
                {
                    newZeroAreaTriangles.push_back(newTriangle_0);
                }
                if (!triangleHasZeroLengthEdges(newTriangle_1, eps)
                    && !containedAndErase(newTriangle_1, newZeroAreaTriangles, eps))
                {
                    newZeroAreaTriangles.push_back(newTriangle_1);
                }
            } else {
                newTriangles.push_back(newTriangle_0);
                newTriangles.push_back(newTriangle_1);
            }
            return 1;
        }
    }

    // Only one vertex belongs to the cut plane.
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) != 3 ? (i + 1) : 0;
        int k = (j + 1) != 3 ? (j + 1) : 0;

        if (vertexTags[i] == 0) {
            if (vertexTags[j] == 1 && vertexTags[k] == 1) {
                if (triangleHasZeroArea(triangle, eps)) {
                    if (!triangleHasZeroLengthEdges(triangle, eps)
                        && !containedAndErase(triangle, newZeroAreaTriangles, eps))
                    {
                        newZeroAreaTriangles.push_back(triangle);
                    }
                } else {
                    newTriangles.push_back(triangle);
                }
                return 1;
            } else if (vertexTags[j] == -1 && vertexTags[k] == -1) {
                return 0;
            } else if (vertexTags[j] == 1 && vertexTags[k] == -1) {
                Array<T, 3> intersection((T)0.0, (T)0.0, (T)0.0);
                int rv = 0;
                rv = lineIntersectionWithPlane<T>(
                    plane, triangle[j], triangle[k], eps, intersection);
                if (rv != 1) {
                    return -1;
                }
                Triangle newTriangle(triangle[i], triangle[j], intersection);
                if (triangleHasZeroArea(triangle, eps)) {
                    if (!triangleHasZeroLengthEdges(newTriangle, eps)
                        && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                    {
                        newZeroAreaTriangles.push_back(newTriangle);
                    }
                } else {
                    newTriangles.push_back(newTriangle);
                }
                return 1;
            } else if (vertexTags[j] == -1 && vertexTags[k] == 1) {
                Array<T, 3> intersection((T)0.0, (T)0.0, (T)0.0);
                int rv = 0;
                rv = lineIntersectionWithPlane<T>(
                    plane, triangle[j], triangle[k], eps, intersection);
                if (rv != 1) {
                    return -1;
                }
                Triangle newTriangle(triangle[i], intersection, triangle[k]);
                if (triangleHasZeroArea(triangle, eps)) {
                    if (!triangleHasZeroLengthEdges(newTriangle, eps)
                        && !containedAndErase(newTriangle, newZeroAreaTriangles, eps))
                    {
                        newZeroAreaTriangles.push_back(newTriangle);
                    }
                } else {
                    newTriangles.push_back(newTriangle);
                }
                return 1;
            }
        }
    }

    // Only two of the three vertices belong to the cut plane.
    for (int i = 0; i < 3; i++) {
        int j = (i + 1) != 3 ? (i + 1) : 0;
        int k = (j + 1) != 3 ? (j + 1) : 0;

        if (vertexTags[i] == 0 && vertexTags[j] == 0) {
            if (vertexTags[k] == 1) {
                if (triangleHasZeroArea(triangle, eps)) {
                    if (!triangleHasZeroLengthEdges(triangle, eps)
                        && !containedAndErase(triangle, newZeroAreaTriangles, eps))
                    {
                        newZeroAreaTriangles.push_back(triangle);
                    }
                } else {
                    newTriangles.push_back(triangle);
                }
                return 1;
            } else if (vertexTags[k] == -1) {
                return 0;
            }
        }
    }

    // All 3 vertices belong to the cut plane.
    if (vertexTags[0] == 0 && vertexTags[1] == 0 && vertexTags[2] == 0) {
        if (triangleHasZeroArea(triangle, eps)) {
            if (!triangleHasZeroLengthEdges(triangle, eps)
                && !containedAndErase(triangle, newZeroAreaTriangles, eps))
            {
                newZeroAreaTriangles.push_back(triangle);
            }
        } else {
            newTriangles.push_back(triangle);
        }
        return 1;
    }

    return -1;
}

template <typename T>
int TriangleSet<T>::cutWithPlane(Plane<T> const &plane, TriangleSet<T> &newTriangleSet) const
{
    T norm_normal = norm(plane.normal);
    PLB_ASSERT(!util::isZero(norm_normal, eps));
    Plane<T> newPlane;
    newPlane.point = plane.point;
    newPlane.normal = plane.normal / norm_normal;

    newTriangleSet.triangles.clear();
    newTriangleSet.eps = eps;
    newTriangleSet.minEdgeLength = std::numeric_limits<T>::max();
    newTriangleSet.maxEdgeLength = std::numeric_limits<T>::min();
    newTriangleSet.minTriangleArea = std::numeric_limits<T>::max();
    newTriangleSet.maxTriangleArea = std::numeric_limits<T>::min();
    T huge = std::numeric_limits<T>::max();
    newTriangleSet.boundingCuboid.lowerLeftCorner = Array<T, 3>((T)huge, (T)huge, (T)huge);
    newTriangleSet.boundingCuboid.upperRightCorner = Array<T, 3>((T)-huge, (T)-huge, (T)-huge);

    std::vector<Triangle> newTriangles;
    std::vector<Triangle> newZeroAreaTriangles;
    for (pluint iTriangle = 0; iTriangle < triangles.size(); iTriangle++) {
        if (cutTriangleWithPlane(newPlane, triangles[iTriangle], newTriangles, newZeroAreaTriangles)
            == -1) {
            return -1;
        }
    }
    newTriangles.insert(
        newTriangles.end(), newZeroAreaTriangles.begin(), newZeroAreaTriangles.end());

    newTriangleSet.triangles.swap(newTriangles);

    if (newTriangleSet.triangles.size() != 0) {
        newTriangleSet.computeMinMaxEdges();
        newTriangleSet.computeMinMaxAreas();
        newTriangleSet.computeBoundingCuboid();
    }

    if (newTriangleSet.triangles.size() == 0 || newTriangleSet.triangles.size() == triangles.size())
    {
        return 0;
    }

    return 1;
}

template <typename T>
int TriangleSet<T>::cutWithPlane(
    Plane<T> const &plane, Cuboid<T> const &cuboid, TriangleSet<T> &newTriangleSet) const
{
    T norm_normal = norm(plane.normal);
    PLB_ASSERT(!util::isZero(norm_normal, eps));
    Plane<T> newPlane;
    newPlane.point = plane.point;
    newPlane.normal = plane.normal / norm_normal;

#ifdef PLB_DEBUG
    T norm_diagonal = norm(cuboid.upperRightCorner - cuboid.lowerLeftCorner);
#endif
    PLB_ASSERT(!util::isZero(norm_diagonal, eps));

    newTriangleSet.triangles.clear();
    newTriangleSet.eps = eps;
    newTriangleSet.minEdgeLength = std::numeric_limits<T>::max();
    newTriangleSet.maxEdgeLength = std::numeric_limits<T>::min();
    newTriangleSet.minTriangleArea = std::numeric_limits<T>::max();
    newTriangleSet.maxTriangleArea = std::numeric_limits<T>::min();
    T huge = std::numeric_limits<T>::max();
    newTriangleSet.boundingCuboid.lowerLeftCorner = Array<T, 3>((T)huge, (T)huge, (T)huge);
    newTriangleSet.boundingCuboid.upperRightCorner = Array<T, 3>((T)-huge, (T)-huge, (T)-huge);

    std::vector<Triangle> newTriangles;
    std::vector<Triangle> newZeroAreaTriangles;

    for (pluint iTriangle = 0; iTriangle < triangles.size(); iTriangle++) {
        Triangle const &triangle = triangles[iTriangle];

        Array<T, 3> vertices[3];
        vertices[0] = triangle[0];
        vertices[1] = triangle[1];
        vertices[2] = triangle[2];

        // Check if the triangle is fully contained in the cuboid.
        int isNotFullyContained = 0;
        for (int iVertex = 0; iVertex < 3; iVertex++) {
            Array<T, 3> diff_l;
            diff_l = vertices[iVertex] - cuboid.lowerLeftCorner;

            Array<T, 3> diff_u;
            diff_u = vertices[iVertex] - cuboid.upperRightCorner;

            if (util::lessThan(diff_l[0], (T)0, eps) || util::lessThan(diff_l[1], (T)0, eps)
                || util::lessThan(diff_l[2], (T)0, eps) || util::greaterThan(diff_u[0], (T)0, eps)
                || util::greaterThan(diff_u[1], (T)0, eps)
                || util::greaterThan(diff_u[2], (T)0, eps))
            {
                isNotFullyContained = 1;
                break;
            }
        }

        if (isNotFullyContained) {
            if (triangleHasZeroArea(triangle, eps)) {
                if (!triangleHasZeroLengthEdges(triangle, eps)
                    && !containedAndErase(triangle, newZeroAreaTriangles, eps))
                {
                    newZeroAreaTriangles.push_back(triangle);
                }
            } else {
                newTriangles.push_back(triangle);
            }
            continue;
        }

        if (cutTriangleWithPlane(newPlane, triangle, newTriangles, newZeroAreaTriangles) == -1)
            return -1;
    }

    newTriangles.insert(
        newTriangles.end(), newZeroAreaTriangles.begin(), newZeroAreaTriangles.end());

    newTriangleSet.triangles.swap(newTriangles);

    if (newTriangleSet.triangles.size() != 0) {
        newTriangleSet.computeMinMaxEdges();
        newTriangleSet.computeMinMaxAreas();
        newTriangleSet.computeBoundingCuboid();
    }

    if (newTriangleSet.triangles.size() == 0 || newTriangleSet.triangles.size() == triangles.size())
    {
        return 0;
    }

    return 1;
}

template <typename T>
Array<T, 3> TriangleSet<T>::getCentroid() const
{
    plint size = triangles.size();
    Array<T, 3> centroid((T)0, (T)0, (T)0);
    plint n = 0;
    for (plint i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            centroid += triangles[i][j];
            n++;
        }
    }

    if (n != 0) {
        centroid /= (T)n;
    }
    return centroid;
}

/// IMPLEMENTED BY HELEN MORRISON (2015)
// Based on
// http://stackoverflow.com/questions/2083771/a-method-to-calculate-the-centre-of-mass-from-a-stl-stereo-lithography-file
// (last accessed: 2015-06-04 18:43)
template <typename T>
Array<T, 3> TriangleSet<T>::getCenterOfMass() const
{
    int numTriangles = triangles.size();

    T totalVolume = 0, currentVolume;
    T xCenter = 0, yCenter = 0, zCenter = 0;

    for (int i = 0; i < numTriangles; i++) {
        currentVolume = (triangles[i][0][0] * triangles[i][1][1] * triangles[i][2][2]
                         - triangles[i][0][0] * triangles[i][2][1] * triangles[i][1][2]
                         - triangles[i][1][0] * triangles[i][0][1] * triangles[i][2][2]
                         + triangles[i][1][0] * triangles[i][2][1] * triangles[i][0][2]
                         + triangles[i][2][0] * triangles[i][0][1] * triangles[i][1][2]
                         - triangles[i][2][0] * triangles[i][1][1] * triangles[i][0][2])
                        / 6;
        totalVolume += currentVolume;
        xCenter +=
            ((triangles[i][0][0] + triangles[i][1][0] + triangles[i][2][0]) / 4) * currentVolume;
        yCenter +=
            ((triangles[i][0][1] + triangles[i][1][1] + triangles[i][2][1]) / 4) * currentVolume;
        zCenter +=
            ((triangles[i][0][2] + triangles[i][1][2] + triangles[i][2][2]) / 4) * currentVolume;
    }

    return Array<T, 3>(xCenter / totalVolume, yCenter / totalVolume, zCenter / totalVolume);
}

template <typename T>
void TriangleSet<T>::computeMinMaxEdges()
{
    minEdgeLength = std::numeric_limits<T>::max();
    maxEdgeLength = std::numeric_limits<T>::min();

    T nextMin, nextMax;
    for (pluint i = 0; i < triangles.size(); ++i) {
        computeMinMaxEdge(i, nextMin, nextMax);
        minEdgeLength = std::min(minEdgeLength, nextMin);
        maxEdgeLength = std::max(maxEdgeLength, nextMax);
    }
}

template <typename T>
void TriangleSet<T>::computeMinMaxEdge(pluint iTriangle, T &minEdge, T &maxEdge) const
{
    PLB_ASSERT(iTriangle < triangles.size());
    Triangle const &triangle = triangles[iTriangle];
    T edge1 = norm(triangle[1] - triangle[0]);
    T edge2 = norm(triangle[2] - triangle[1]);
    T edge3 = norm(triangle[0] - triangle[2]);
    minEdge = std::min(edge1, std::min(edge2, edge3));
    maxEdge = std::max(edge1, std::max(edge2, edge3));
}

template <typename T>
void TriangleSet<T>::computeMinMaxAreas()
{
    minTriangleArea = std::numeric_limits<T>::max();
    maxTriangleArea = std::numeric_limits<T>::min();

    for (pluint i = 0; i < triangles.size(); ++i) {
        T area = computeArea(i);
        minTriangleArea = std::min(minTriangleArea, area);
        maxTriangleArea = std::max(maxTriangleArea, area);
    }
}

template <typename T>
T TriangleSet<T>::computeArea(pluint iTriangle) const
{
    PLB_ASSERT(iTriangle < triangles.size());
    return computeTriangleArea(
        triangles[iTriangle][0], triangles[iTriangle][1], triangles[iTriangle][2]);
}

template <typename T>
void TriangleSet<T>::computeBoundingCuboid()
{
    T xMin, yMin, zMin;
    T xMax, yMax, zMax;

    xMin = yMin = zMin = std::numeric_limits<T>::max();
    xMax = yMax = zMax = -std::numeric_limits<T>::max();
    for (pluint i = 0; i < triangles.size(); ++i) {
        Triangle const &triangle = triangles[i];

        xMin = std::min(xMin, std::min(triangle[0][0], std::min(triangle[1][0], triangle[2][0])));
        yMin = std::min(yMin, std::min(triangle[0][1], std::min(triangle[1][1], triangle[2][1])));
        zMin = std::min(zMin, std::min(triangle[0][2], std::min(triangle[1][2], triangle[2][2])));

        xMax = std::max(xMax, std::max(triangle[0][0], std::max(triangle[1][0], triangle[2][0])));
        yMax = std::max(yMax, std::max(triangle[0][1], std::max(triangle[1][1], triangle[2][1])));
        zMax = std::max(zMax, std::max(triangle[0][2], std::max(triangle[1][2], triangle[2][2])));
    }
    boundingCuboid.lowerLeftCorner = Array<T, 3>(xMin, yMin, zMin);
    boundingCuboid.upperRightCorner = Array<T, 3>(xMax, yMax, zMax);
}

template <typename T>
void TriangleSet<T>::skipLines(plint nLines, FILE *fp) const
{
    for (plint i = 0; i < nLines; i++)
        while (fgetc(fp) != '\n')
            ;
}

template <typename T>
int TriangleSet<T>::readAhead(FILE *fp, char commentCharacter) const
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
#else
            (void)ungetc(nextChar, fp);
#endif
            PLB_ASSERT(rv != EOF);  // Unexpected error.
            return 0;
        }
    }
    return EOF;
}

template <typename T>
bool TriangleSet<T>::checkForBufferOverflow(char *buf) const
{
    while (*buf != '\0') {
        if (*buf == '\n') {
            return true;  // All the line was read into the buffer.
        }
        buf++;
    }
    return false;  // The full line of text was not read into the buffer.
}

}  // namespace plb

#undef PLB_CBUFSIZ

#endif  // TRIANGLE_SET_HH
