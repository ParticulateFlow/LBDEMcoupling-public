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

#ifndef META_STUFF_HELPER_H
#define META_STUFF_HELPER_H

#include <map>
#include <set>

#include "atomicBlock/atomicContainerBlock2D.h"

namespace plb {

struct VectorIsLess {
    bool operator()(std::vector<int> const &v1, std::vector<int> const &v2) const
    {
        pluint bound = std::max(v1.size(), v2.size());
        for (pluint i = 0; i < bound; ++i) {
            int val1 = i < v1.size() ? v1[i] : -1;
            int val2 = i < v2.size() ? v2[i] : -1;
            if (val1 < val2) {
                return true;
            } else if (val1 > val2) {
                return false;
            }
        }
        return false;
    }
};

inline bool vectorEquals(std::vector<int> const &v1, std::vector<int> const &v2)
{
    return !(VectorIsLess()(v1, v2) || VectorIsLess()(v2, v1));
}

/// Container object for StoreDynamicsFunctionalXD
class StoreDynamicsID : public ContainerBlockData {
public:
    typedef std::set<std::vector<int>, VectorIsLess> ChainCollection;

public:
    virtual StoreDynamicsID *clone() const
    {
        return new StoreDynamicsID(*this);
    }
    void addIdChain(std::vector<int> const &chain)
    {
        idChains.insert(chain);
    }
    ChainCollection const &getIds() const
    {
        return idChains;
    }
    void startIterations()
    {
        pos = idChains.rbegin();
    }
    std::vector<int> iterate()
    {
        ++pos;
        if (empty()) {
            std::vector<int> none;
            none.push_back(-1);
            return none;
        } else {
            return *pos;
        }
    }
    bool empty() const
    {
        return pos == idChains.rend();
    }
    std::vector<int> const &getCurrent() const
    {
        return *pos;
    }

private:
    ChainCollection::const_reverse_iterator pos;
    ChainCollection idChains;
};

}  // namespace plb

#endif  // META_STUFF_HELPER_H
