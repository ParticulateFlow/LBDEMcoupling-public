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

#ifndef CO_PROCESSOR_COMMUNICATION_3D_HH
#define CO_PROCESSOR_COMMUNICATION_3D_HH

#include <iomanip>

#include "coProcessors/coProcessorCommunication3D.h"
#include "core/globalDefs.h"
#include "io/parallelIO.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
void transferToCoProcessors(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    MultiBlockManagement3D const &management = lattice.getMultiBlockManagement();
    ThreadAttribution const &threadAttribution = management.getThreadAttribution();

    std::vector<char> data;
    for (pluint iBlock = 0; iBlock < management.getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = management.getLocalInfo().getBlocks()[iBlock];
        plint handle = threadAttribution.getCoProcessorHandle(blockId);
        if (handle >= 0) {
            BlockLattice3D<T, Descriptor> &component = lattice.getComponent(blockId);
            component.getDataTransfer().send(
                component.getBoundingBox(), data, modif::staticVariables);
            global::defaultCoProcessor3D<T>().send(handle, component.getBoundingBox(), data);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void transferFromCoProcessors(MultiBlockLattice3D<T, Descriptor> &lattice)
{
    MultiBlockManagement3D const &management = lattice.getMultiBlockManagement();
    ThreadAttribution const &threadAttribution = management.getThreadAttribution();

    std::vector<char> data;
    for (pluint iBlock = 0; iBlock < management.getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = management.getLocalInfo().getBlocks()[iBlock];
        plint handle = threadAttribution.getCoProcessorHandle(blockId);
        if (handle >= 0) {
            BlockLattice3D<T, Descriptor> &component = lattice.getComponent(blockId);
            global::defaultCoProcessor3D<T>().receive(handle, component.getBoundingBox(), data);
            component.getDataTransfer().receive(
                component.getBoundingBox(), data, modif::staticVariables);
        }
    }
}

template <typename T, template <typename U> class Descriptor>
void initiateCoProcessors(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint dynamicsId, bool printInfo)
{
    std::map<plint, PureDynamics<T> > dynamicsPattern =
        identifyBlocksWithPureDynamics(lattice, dynamicsId);

    if (printInfo) {
        typename std::map<plint, PureDynamics<T> >::const_iterator it = dynamicsPattern.begin();
        plint i = 0;
        for (; it != dynamicsPattern.end(); ++it, ++i) {
            plint blockId = it->first;
            bool isPure = it->second.isPure;
            T omega = it->second.omega;
            if (isPure) {
                pcout << std::setw(4) << blockId << ": (" << std::setw(8) << omega << "); ";
            } else {
                pcout << std::setw(4) << blockId << ";       ";
            }
            if (i % 6 == 0)
                pcout << std::endl;
        }
        pcout << std::endl;
    }

    MultiBlockManagement3D const &management = lattice.getMultiBlockManagement();
    std::map<plint, int> coprocessors;
    for (pluint iBlock = 0; iBlock < management.getLocalInfo().getBlocks().size(); ++iBlock) {
        plint blockId = management.getLocalInfo().getBlocks()[iBlock];

        typename std::map<plint, PureDynamics<T> >::const_iterator itPure =
            dynamicsPattern.find(blockId);
        PLB_ASSERT(itPure != dynamicsPattern.end());
        bool isPure = itPure->second.isPure;
        T omega = itPure->second.omega;
        if (isPure) {
            SmartBulk3D smartBulk(management, blockId);
            Box3D domain(smartBulk.computeEnvelope());
            int handle;
            global::defaultCoProcessor3D<T>().addDomain(
                domain.getNx(), domain.getNy(), domain.getNz(), omega, handle);
            coprocessors.insert(std::pair<plint, int>(blockId, handle));
        } else {
            coprocessors.insert(std::pair<plint, int>(blockId, -1));
        }
    }
    lattice.setCoProcessors(coprocessors);
}

}  // namespace plb

#endif  // CO_PROCESSOR_COMMUNICATION_3D_HH
