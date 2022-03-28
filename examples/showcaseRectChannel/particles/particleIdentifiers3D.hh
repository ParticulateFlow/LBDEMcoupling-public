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

#ifndef PARTICLE_IDENTIFIERS_3D_HH
#define PARTICLE_IDENTIFIERS_3D_HH

#include <sstream>

#include "core/blockIdentifiers.h"
#include "core/runTimeDiagnostics.h"

namespace plb {

namespace meta {

template <typename T, template <typename U> class Descriptor>
ParticleRegistration3D<T, Descriptor>::~ParticleRegistration3D()
{
    for (pluint iEntry = 0; iEntry < particleByNumber.size(); ++iEntry) {
        delete particleByNumber[iEntry].generator;
    }
}

template <typename T, template <typename U> class Descriptor>
int ParticleRegistration3D<T, Descriptor>::announce(
    std::string nameOfParticle, ParticleGenerator3D<T, Descriptor> *generator)
{
    Entry entry(nameOfParticle, generator);
    typename EntryMap::iterator it = particleByName.find(entry);
    if (it != particleByName.end()) {
        plbLogicError(
            std::string("The particle class ") + nameOfParticle
            + std::string(" was registered twice"));
    }
    particleByNumber.push_back(entry);
    int nextId = particleByNumber.size();
    particleByName[entry] = nextId;
    return nextId;
}

template <typename T, template <typename U> class Descriptor>
int ParticleRegistration3D<T, Descriptor>::getId(std::string name) const
{
    Entry entry(name, 0);
    typename EntryMap::const_iterator it = particleByName.find(entry);
    if (it == particleByName.end()) {
        return 0;
    } else {
        return it->second;
    }
}

template <typename T, template <typename U> class Descriptor>
int ParticleRegistration3D<T, Descriptor>::getNumId() const
{
    return (int)(particleByNumber.size());
}

template <typename T, template <typename U> class Descriptor>
std::string ParticleRegistration3D<T, Descriptor>::getName(int id) const
{
    if (id == 0) {
        return std::string("Undefined");
    }
    if (id < 0 || id > (int)particleByNumber.size()) {
        std::stringstream message;
        message << "A particle class with ID " << id << " doesn't exist.";
        plbLogicError(message.str());
    }
    return particleByNumber[id - 1].name;
}

template <typename T, template <typename U> class Descriptor>
Particle3D<T, Descriptor> *ParticleRegistration3D<T, Descriptor>::generate(
    HierarchicUnserializer &unserializer)
{
    plint id = unserializer.getId();
    PLB_ASSERT(id > 0 && (pluint)id <= particleByNumber.size());
    return particleByNumber[id - 1].generator->generate(unserializer);
}

template <typename T, template <typename U> class Descriptor>
typename ParticleRegistration3D<T, Descriptor>::EntryMap::const_iterator
    ParticleRegistration3D<T, Descriptor>::begin() const
{
    return particleByName.begin();
}

template <typename T, template <typename U> class Descriptor>
typename ParticleRegistration3D<T, Descriptor>::EntryMap::const_iterator
    ParticleRegistration3D<T, Descriptor>::end() const
{
    return particleByName.end();
}

template <typename T, template <typename U> class Descriptor>
ParticleRegistration3D<T, Descriptor> &particleRegistration3D()
{
    static ParticleRegistration3D<T, Descriptor> instance;
    return instance;
}

}  // namespace meta

}  // namespace plb

#endif  // PARTICLE_IDENTIFIERS_3D_HH
