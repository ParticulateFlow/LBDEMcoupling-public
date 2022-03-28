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

#ifndef STL_FILE_IO_H
#define STL_FILE_IO_H

/*
 * STL I/O stuff. Everything is independent of whether it's going to be used
 * on a RawTriangleMesh, or a RawConnectedTriangleMesh, or a DEFtriangleMesh.
 */

#include <string>
#include <vector>

#include "core/globalDefs.h"
#include "offLattice/triangleMesh.h"
#include "offLattice/triangleSelector.h"

namespace plb {

template <typename T>
class STLreader {
public:
    typedef Array<Array<T, 3>, 3> Triangle;
    STLreader(std::string fname, T eps_ = getEpsilon<T>(DBL), TriangleSelector<T> *selector = 0);
    std::vector<std::vector<Triangle> > const &getParts() const
    {
        return parts;
    }
    std::vector<std::string> const &getNames() const
    {
        return names;
    }
    T getEps() const
    {
        return eps;
    }

private:
    void readSTL(std::string fname, TriangleSelector<T> *selector);
    bool isAsciiSTL(FILE *fp, std::string fname);
    bool readAsciiSTL(FILE *fp, TriangleSelector<T> *selector);
    bool readBinarySTL(FILE *fp, TriangleSelector<T> *selector);
    bool triangleHasZeroLengthEdges(Triangle const &triangle) const;
    void fixOrientation(Triangle &triangle, Array<T, 3> const &n) const;
    bool checkForBufferOverflow(char *buf);
    char *cleanString(char *str);

private:
    std::vector<std::vector<Triangle> > parts;
    std::vector<std::string> names;
    T eps;
};

// In the following functions, if mainProcOnly = false, then all processes write to disk. This fact
// must be taken under consideration by the caller when the fname parameter is created.

template <typename T>
void writeAsciiSTL(
    TriangleMesh<T> &mesh, std::string fname, int numDecimalDigits = 10, bool mainProcOnly = true);

template <typename T>
void writeMultiPartAsciiSTL(
    TriangleMesh<T> &mesh, std::string fname, int numDecimalDigits = 10, bool mainProcOnly = true);

template <typename T>
void writeBinarySTL(
    typename TriangleMesh<T>::PTriangleIterator const &iterator, std::string fname,
    bool mainProcOnly = true);

template <typename T>
void writeBinarySTL(TriangleMesh<T> &mesh, std::string fname, bool mainProcOnly = true);

template <typename T>
void writeMultiPartBinarySTL(TriangleMesh<T> &mesh, std::string fname, bool mainProcOnly = true);

}  // namespace plb

#endif  // STL_FILE_IO_H
