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

#include "io/vtkStructuredDataOutput.h"

#include "core/globalDefs.h"
#include "core/serializer.h"
#include "core/serializer.hh"
#include "io/base64.h"
#include "io/base64.hh"
#include "io/serializerIO.h"
#include "io/serializerIO.hh"
#include "io/vtkStructuredDataOutput.hh"
#include "parallelism/mpiManager.h"

namespace plb {

////////// class VtkStructuredWriter3D ////////////////////////////////////////

VtkStructuredWriter3D::VtkStructuredWriter3D(std::string const &fileName_) :
    fileName(fileName_), ostr(0)
{
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fileName.c_str());
        if (!(*ostr)) {
            std::cerr << "could not open file " << fileName << "\n";
            return;
        }
    }
}

VtkStructuredWriter3D::~VtkStructuredWriter3D()
{
    delete ostr;
}

void VtkStructuredWriter3D::writeHeader(Box3D domain)
{
    if (global::mpi().isMainProcessor()) {
        (*ostr) << "<?xml version=\"1.0\"?>\n";
#ifdef PLB_BIG_ENDIAN
        (*ostr) << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
#else
        (*ostr)
            << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
#endif
        (*ostr) << "<StructuredGrid WholeExtent=\"" << domain.x0 << " " << domain.x1 << " "
                << domain.y0 << " " << domain.y1 << " " << domain.z0 << " " << domain.z1 << "\">\n";
    }
}

void VtkStructuredWriter3D::startPiece(Box3D domain, const Array<double, 3> &origin, double deltaX)
{
    if (global::mpi().isMainProcessor()) {
        (*ostr) << "<Piece Extent=\"" << domain.x0 << " " << domain.x1 << " " << domain.y0 << " "
                << domain.y1 << " " << domain.z0 << " " << domain.z1 << "\">\n";
        (*ostr) << "<Points>\n";
        (*ostr) << "<DataArray  NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">\n";
        // point loop
        plb::plint i, j, k;
        for (i = domain.z0; i <= domain.z1; i++) {
            for (j = domain.y0; j <= domain.y1; j++) {
                for (k = domain.x0; k <= domain.x1; k++) {
                    (*ostr) << k * deltaX - origin[2] << " " << j * deltaX - origin[1] << " "
                            << i * deltaX - origin[0] << " \n";
                }
            }
        }
        (*ostr) << "</DataArray>\n</Points>\n";
        (*ostr) << "<PointData>\n";
    }
}

void VtkStructuredWriter3D::endPiece()
{
    if (global::mpi().isMainProcessor()) {
        (*ostr) << "</PointData>\n";
        (*ostr) << "</Piece>\n";
    }
}

void VtkStructuredWriter3D::writeFooter()
{
    if (global::mpi().isMainProcessor()) {
        (*ostr) << "</StructuredGrid>\n";
        (*ostr) << "</VTKFile>\n";
    }
}

}  // namespace plb
