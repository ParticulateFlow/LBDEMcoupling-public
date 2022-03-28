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

#ifndef CO_PROCESSOR_INSTANTIATION_3D_HH
#define CO_PROCESSOR_INSTANTIATION_3D_HH

#include "coProcessors/coProcessorInstantiation3D.h"
#include "core/globalDefs.h"

namespace plb {

template <typename T, template <typename U> class Descriptor>
std::map<plint, PureDynamics<T> > identifyBlocksWithPureDynamics(
    MultiBlockLattice3D<T, Descriptor> &lattice, plint dynamicsId)
{
    MultiContainerBlock3D container(lattice);
    std::vector<MultiBlock3D *> latticeAndContainer;
    latticeAndContainer.push_back(&lattice);
    latticeAndContainer.push_back(&container);

    applyProcessingFunctional(
        new IdentifyPureDynamics3D<T, Descriptor>(dynamicsId), lattice.getBoundingBox(),
        latticeAndContainer);

    MultiBlockManagement3D const &management = container.getMultiBlockManagement();
    SparseBlockStructure3D const &sparseBlock = management.getSparseBlockStructure();
    ThreadAttribution const &threadAttribution = management.getThreadAttribution();
    std::map<plint, Box3D> const &domains = sparseBlock.getBulks();

    std::vector<plint> domainIds(domains.size());
    std::vector<int> isPure(domains.size());
    std::vector<T> omega(domains.size());

    std::map<plint, Box3D>::const_iterator it = domains.begin();
    plint pos = 0;
    for (; it != domains.end(); ++it) {
        plint id = it->first;
        domainIds[pos] = id;
        if (threadAttribution.isLocal(id)) {
            AtomicContainerBlock3D const &component = container.getComponent(id);
            PureDynamics<T> const *data =
                dynamic_cast<PureDynamics<T> const *>(component.getData());
            PLB_ASSERT(data);
            isPure[pos] = data->isPure ? 1 : 0;
            omega[pos] = data->omega;
        } else {
            isPure[pos] = 0;
            omega[pos] = T();
        }
        ++pos;
    }

#ifdef PLB_MPI_PARALLEL
    std::vector<int> isPureTmp(isPure.size());
    global::mpi().reduceVect(isPure, isPureTmp, MPI_SUM);
    global::mpi().bCast(&isPureTmp[0], isPureTmp.size());
    isPureTmp.swap(isPure);

    std::vector<T> omegaTmp(omega.size());
    global::mpi().reduceVect(omega, omegaTmp, MPI_SUM);
    global::mpi().bCast(&omegaTmp[0], omegaTmp.size());
    omegaTmp.swap(omega);
#endif

    std::map<plint, PureDynamics<T> > pureDynamicsMap;
    for (pluint i = 0; i < isPure.size(); ++i) {
        pureDynamicsMap.insert(
            std::pair<plint, PureDynamics<T> >(domainIds[i], PureDynamics<T>(isPure[i], omega[i])));
    }

    return pureDynamicsMap;
}

}  // namespace plb

#endif  // CO_PROCESSOR_INSTANTIATION_3D_HH
