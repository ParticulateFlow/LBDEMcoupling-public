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

#ifndef CO_PROCESSOR_FUNCTIONAL_3D_HH
#define CO_PROCESSOR_FUNCTIONAL_3D_HH

#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "coProcessors/coProcessorFunctional3D.h"
#include "core/plbDebug.h"
#include "core/util.h"

namespace plb {

/* ********  IdentifyPureDynamics3D ************************************ */

template <typename T, template <typename U> class Descriptor>
IdentifyPureDynamics3D<T, Descriptor>::IdentifyPureDynamics3D(plint dynamicsId_) :
    dynamicsId(dynamicsId_)
{ }

template <typename T, template <typename U> class Descriptor>
void IdentifyPureDynamics3D<T, Descriptor>::processGenericBlocks(
    Box3D domain, std::vector<AtomicBlock3D *> blocks)
{
    T eps = 1.e-6;
    PLB_PRECONDITION(blocks.size() == 2);
    BlockLattice3D<T, Descriptor> &lattice =
        *dynamic_cast<BlockLattice3D<T, Descriptor> *>(blocks[0]);
    AtomicContainerBlock3D &container = *dynamic_cast<AtomicContainerBlock3D *>(blocks[1]);
    PureDynamics<T> *pureDynamics = new PureDynamics<T>;
    bool firstElement = true;
    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                Dynamics<T, Descriptor> &dynamics = lattice.get(iX, iY, iZ).getDynamics();
                pureDynamics->isPure = dynamics.getId() == dynamicsId;
                if (firstElement) {
                    pureDynamics->omega = dynamics.getOmega();
                } else {
                    if (std::fabs(pureDynamics->omega - dynamics.getOmega()) > eps) {
                        pureDynamics->isPure = false;
                    }
                }
                if (!pureDynamics->isPure)
                    break;
                firstElement = false;
            }
            if (!pureDynamics->isPure)
                break;
        }
        if (!pureDynamics->isPure)
            break;
    }
    if (!pureDynamics->isPure) {
        pureDynamics->omega = T();
    }
    container.setData(pureDynamics);
}

template <typename T, template <typename U> class Descriptor>
IdentifyPureDynamics3D<T, Descriptor> *IdentifyPureDynamics3D<T, Descriptor>::clone() const
{
    return new IdentifyPureDynamics3D<T, Descriptor>(*this);
}

template <typename T, template <typename U> class Descriptor>
void IdentifyPureDynamics3D<T, Descriptor>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template <typename T, template <typename U> class Descriptor>
BlockDomain::DomainT IdentifyPureDynamics3D<T, Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // CO_PROCESSOR_FUNCTIONAL_3D_HH
