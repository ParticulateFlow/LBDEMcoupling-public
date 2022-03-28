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

#ifdef HDF5

#ifndef XDMF_DATA_OUTPUT_H
#define XDMF_DATA_OUTPUT_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataField3D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "core/serializer.h"
#include "io/hdfWrapper.h"
#include "io/multiBlockWriter3D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

class ParallelXdmfDataWriter3D {
private:
    std::string xdmf_fname;
    std::string h5_fname;
    std::ofstream *fhandle;
    int field = 0;

public:
    ParallelXdmfDataWriter3D(std::string fname) :
        xdmf_fname(FileName(fname + ".xdmf").defaultPath(global::directories().getOutputDir())),
        h5_fname(FileName(fname + ".h5").defaultPath(global::directories().getOutputDir()))
    {
        if (global::mpi().isMainProcessor()) {
            fhandle = new std::ofstream(xdmf_fname.c_str());
            if (!(*fhandle)) {
                std::cerr << "could not open file " << fname << "\n";
                return;
            }
            (*fhandle) << "<?xml version = \"1.0\"?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
            (*fhandle) << "\t<Xdmf Version=\"2.0\">\n";
            (*fhandle) << "\t\t<Domain>\n";
            (*fhandle).flush();
        }
    }
    ~ParallelXdmfDataWriter3D()
    {
        if (global::mpi().isMainProcessor()) {
            (*fhandle) << "\t\t</Domain>\n";
            (*fhandle) << "\t</Xdmf>\n";
            (*fhandle).close();
            delete fhandle;
        }
    }
    template <typename T>
    void writeDataField(MultiBlock3D &multiBlock, const std::string &field_name)
    {
        std::vector<plint> offset;
        std::vector<plint> myBlockIds;
        std::vector<std::vector<char>> data;

        MultiBlockManagement3D const &management = multiBlock.getMultiBlockManagement();
        std::map<plint, Box3D> const &bulks = management.getSparseBlockStructure().getBulks();

        parallelIO::dumpData(multiBlock, true, offset, myBlockIds, data);
        int total_size = offset.back() / sizeof(T);
        char set_name[20];
        sprintf(set_name, "binary_blob%d", field);

        if (global::mpi().isMainProcessor()) {
            plint pos = 0;
            std::map<plint, Box3D>::const_iterator it = bulks.begin();

            for (; it != bulks.end(); ++it) {
                plint blockId = it->first;
                SmartBulk3D bulk(management, blockId);
                Box3D localBulk(bulk.toLocal(bulk.getBulk()));
                // localBulk is the real domain without the gho
                int dx = it->second.x0;
                int dy = it->second.y0;
                int dz = it->second.z0;

                int start = 0;
                int size = 0;

                if (pos > 0) {
                    size = (offset[pos] - offset[pos - 1]) / sizeof(T);
                    start = offset[pos - 1] / sizeof(T);
                } else {
                    size = offset[pos] / sizeof(T);
                    start = 0;
                }
                int xn = localBulk.getNx(), yn = localBulk.getNy(), zn = localBulk.getNz();
                int scalar_size = xn * yn * zn;
                int dimension = size / scalar_size;
                std::string attribute_type("Scalar");
                if (dimension > 1)
                    attribute_type = "Vector";

                (*fhandle) << "\t\t<Grid Name = \"Grid" << pos << "\">\n";
                (*fhandle) << "\t\t\t<Topology Dimensions = \"" << xn + 1 << " " << yn + 1 << " "
                           << zn + 1 << "\" Type = \"3DCoRectMesh\" />\n";
                (*fhandle) << "\t\t\t<Geometry Origin = \"\" Type = \"ORIGIN_DXDYDZ\">\n";
                (*fhandle) << "\t\t\t\t<DataItem DataType = \"Float\" Dimensions = \"3\" Format = "
                              "\"XML\" >"
                           << dx << " " << dy << " " << dz << "</DataItem>\n";
                (*fhandle) << "\t\t\t\t<DataItem DataType = \"Float\" Dimensions = \"3\" Format = "
                              "\"XML\" >1 1 1</DataItem>\n";
                (*fhandle) << "\t\t\t</Geometry>\n";
                (*fhandle) << "\t\t\t<Attribute Center = \"Cell\" Active = \"1\"  Name = \""
                           << field_name << "\" Type = \"" << attribute_type << "\" >\n";
                (*fhandle) << "\t\t\t\t<DataItem ItemType = \"HyperSlab\" Dimensions = \"" << xn
                           << " " << yn << " " << zn << " " << dimension
                           << "\" Type = \"HyperSlab\" Format = \"XML\"  DataType = \"Float\"  "
                              "Precision = \""
                           << sizeof(T) << "\">\n";
                (*fhandle) << "\t\t\t\t\t<DataItem  Dimensions = \"3\" Format = \"XML\">" << start
                           << "  1  " << size << "</DataItem>\n";
                (*fhandle) << "\t\t\t\t\t<DataItem DataType = \"Float\" ItemType = \"Uniform\" "
                              "Precision = \""
                           << sizeof(T) << "\" Dimensions = \"" << total_size
                           << "\" Format = \"HDF\">" << h5_fname << ":/" << set_name
                           << "</DataItem>\n";
                (*fhandle) << "\t\t\t\t</DataItem>\n";
                (*fhandle) << "\t\t\t</Attribute>\n";
                (*fhandle) << "\t\t</Grid>\n";

                pos++;
            }
            (*fhandle).flush();
        }

        writeParallelHDF5<T>(
            h5_fname.c_str(), set_name, myBlockIds, offset, data, "", global::mpi().getRank(),
            global::mpi().getGlobalCommunicator(), field > 0);
        field++;
    }
};

}  // namespace plb

#endif  // XDMF_DATA_OUTPUT_H

#endif
