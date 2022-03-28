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

#ifndef OFF_FILE_IO_H
#define OFF_FILE_IO_H

/*
 * OFF I/O stuff. Everything is independent of whether it's going to be used
 * on a RawTriangleMesh, or a RawConnectedTriangleMesh, or a DEFtriangleMesh.
 */

#include <vector>

#include "core/array.h"
#include "core/globalDefs.h"
#include "offLattice/triangleMesh.h"

namespace plb {

template <typename T>
class OFFreader {
public:
    OFFreader(std::string fname);
    std::vector<Array<T, 3> > const &getVertices() const
    {
        return vertices;
    }
    std::vector<std::vector<plint> > const &getFacets() const
    {
        return facets;
    }

private:
    void readOFF(std::string fname);
    // Caution: in the readAsciiOFF method we keep all information that resides
    //          in the OFF file. In the readAsciiOFF method of the TriangleSet class
    //          we neglect all triangles that have one or more edges with length
    //          equal to 0. We do not do this here, because from the OFFreader we
    //          can directly create a connected mesh and we do not want to break
    //          vertex and triangle numbering.
    bool readAsciiOFF(FILE *fp);
    void skipLines(plint nLines, FILE *fp) const;
    int readAhead(FILE *fp, char commentCharacter) const;

private:
    std::vector<Array<T, 3> > vertices;
    std::vector<std::vector<plint> > facets;
};

template <typename T>
void writeAsciiOFF(
    TriangleMesh<T> &mesh, std::string fname, T eps = getEpsilon<T>(DBL),
    int numDecimalDigits = 10);

}  // namespace plb

#endif  // OFF_FILE_IO_H
