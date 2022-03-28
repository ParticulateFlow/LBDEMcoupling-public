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

/** \file
 * Attribution of block id to thread for parallel execution -- header file.
 */

#ifndef THREAD_ATTRIBUTION_H
#define THREAD_ATTRIBUTION_H

#include <map>
#include <vector>

#include "core/globalDefs.h"

namespace plb {

struct ThreadAttribution {
    virtual ~ThreadAttribution() { }
    /// Specifies if the block "blockId" is located on the current MPI thread.
    virtual bool isLocal(plint blockId) const = 0;
    /// True means that all blocks are local, false means we don't know if
    ///   all of them are local or not.
    virtual bool allBlocksAreLocal() const = 0;
    /// Specifies on which MPI process a given block is located.
    virtual int getMpiProcess(plint blockId) const = 0;
    /// Specifies to which of the local shared-memory threads the blockId
    ///   belongs.
    virtual int getLocalThreadId(plint blockId) const = 0;
    /// Merge current attribution with rhs. From rhs, select explicitly
    ///   the specified blocks, and remap them to fit into the new
    ///   attribution, as specified by the parameter remappedIds.
    virtual ThreadAttribution *merge(
        ThreadAttribution const &rhs,
        std::map<plint, std::vector<plint> > const &remappedIds) const = 0;
    /// Extend the current attribution by a selection of new identifiers,
    ///   with explicitly specified mpiProcess and thread attribution.
    virtual ThreadAttribution *extend(
        std::vector<plint> const &ids, std::vector<plint> const &mpiProcesses,
        std::vector<plint> const &localThreads) const = 0;
    virtual bool hasCoProcessors() const
    {
        return false;
    }
    virtual void setCoProcessors(std::map<plint, int> processorHandles) { }
    virtual int getCoProcessorHandle(plint id) const
    {
        return 0;
    }
    virtual ThreadAttribution *clone() const = 0;
};

struct SerialThreadAttribution : public ThreadAttribution {
    virtual bool isLocal(plint blockId) const;
    virtual bool allBlocksAreLocal() const;
    virtual int getMpiProcess(plint blockId) const;
    virtual int getLocalThreadId(plint blockId) const;
    virtual ThreadAttribution *merge(
        ThreadAttribution const &rhs,
        std::map<plint, std::vector<plint> > const &remappedIds) const;
    virtual ThreadAttribution *extend(
        std::vector<plint> const &ids, std::vector<plint> const &mpiProcesses,
        std::vector<plint> const &localThreads) const;
    virtual ThreadAttribution *clone() const;
};

struct OneToOneThreadAttribution : public ThreadAttribution {
    virtual bool isLocal(plint blockId) const;
    virtual bool allBlocksAreLocal() const;
    virtual int getMpiProcess(plint blockId) const;
    virtual int getLocalThreadId(plint blockId) const;
    virtual ThreadAttribution *merge(
        ThreadAttribution const &rhs,
        std::map<plint, std::vector<plint> > const &remappedIds) const;
    virtual ThreadAttribution *extend(
        std::vector<plint> const &ids, std::vector<plint> const &mpiProcesses,
        std::vector<plint> const &localThreads) const;
    virtual ThreadAttribution *clone() const;
};

class ExplicitThreadAttribution : public ThreadAttribution {
public:
    ExplicitThreadAttribution();
    ExplicitThreadAttribution(std::map<plint, plint> const &mpiProcessAttribution_);
    ExplicitThreadAttribution(
        std::map<plint, plint> const &mpiProcessAttribution_,
        std::map<plint, plint> const &localThreadAttribution_);
    void addBlock(plint blockId, plint mpiProcess, plint localThread = 0);
    virtual bool isLocal(plint blockId) const;
    virtual bool allBlocksAreLocal() const;
    virtual int getMpiProcess(plint blockId) const;
    virtual int getLocalThreadId(plint blockId) const;
    virtual ThreadAttribution *merge(
        ThreadAttribution const &rhs,
        std::map<plint, std::vector<plint> > const &remappedIds) const;
    virtual ThreadAttribution *extend(
        std::vector<plint> const &ids, std::vector<plint> const &mpiProcesses,
        std::vector<plint> const &localThreads) const;
    virtual bool hasCoProcessors() const
    {
        return !coProcessors.empty();
    }
    virtual void setCoProcessors(std::map<plint, int> processorHandles)
    {
        coProcessors = processorHandles;
    }
    virtual int getCoProcessorHandle(plint id) const
    {
        std::map<plint, int>::const_iterator it = coProcessors.find(id);
        if (it == coProcessors.end()) {
            return -1;
        } else {
            return it->second;
        }
    }
    virtual ThreadAttribution *clone() const;

private:
    std::map<plint, plint> mpiProcessAttribution;
    std::map<plint, plint> localThreadAttribution;
    std::map<plint, int> coProcessors;
};

}  // namespace plb

#endif  // THREAD_ATTRIBUTION_H
