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

#include "io/serializerIO_3D.h"

#include <fstream>
#include <istream>
#include <ostream>

#include "core/plbDebug.h"
#include "core/runTimeDiagnostics.h"
#include "io/serializerIO.h"
#include "multiBlock/multiBlockSerializer3D.h"
#include "parallelism/mpiManager.h"

namespace plb {

void saveBinaryBlock(Block3D const &block, std::string fName, bool enforceUint)
{
    std::ofstream *ostr = 0;
    bool isOK = true;
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fName.c_str());
        isOK = (bool)(*ostr);
    }
    plbMainProcIOError(
        !isOK, std::string("Could not open binary file ") + fName + std::string(" for saving"));
    serializerToBase64Stream(
        block.getBlockSerializer(
            block.getBoundingBox(), global::IOpolicy().getIndexOrderingForStreams()),
        ostr, enforceUint);
    delete ostr;
}

void loadBinaryBlock(Block3D &block, std::string fName, bool enforceUint)
{
    std::ifstream *istr = 0;
    bool isOK = true;
    if (global::mpi().isMainProcessor()) {
        istr = new std::ifstream(fName.c_str());
        isOK = (bool)(*istr);
    }
    plbMainProcIOError(
        !isOK, std::string("Could not open binary file ") + fName + std::string(" for reading"));
    base64StreamToUnSerializer(
        istr,
        block.getBlockUnSerializer(
            block.getBoundingBox(), global::IOpolicy().getIndexOrderingForStreams()),
        enforceUint);
    delete istr;
}

}  // namespace plb
