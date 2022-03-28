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
 * Geometric operations on collections of 3D domains -- implementation file.
 */

#include "multiBlock/domainManipulation3D.h"

#include "core/plbDebug.h"

namespace plb {

std::vector<DomainAndId3D> getNonOverlapingBlocks(std::vector<DomainAndId3D> const &domainsWithId)
{
    std::vector<DomainAndId3D> nonOverlapingBlocks;
    // Start with the first domain, which is taken without modification.
    if (!domainsWithId.empty()) {
        nonOverlapingBlocks.push_back(domainsWithId[0]);
    }
    // All subsequent domains get special treatment, as their overlap
    //   with previously adopted domains are cut out.
    for (pluint iDomain = 1; iDomain < domainsWithId.size(); ++iDomain) {
        std::vector<Box3D> newDomains;
        newDomains.push_back(domainsWithId[iDomain].domain);
        for (pluint iPrevious = 0; iPrevious < iDomain; ++iPrevious) {
            std::vector<Box3D> exceptedDomains;
            for (pluint iNewPart = 0; iNewPart < newDomains.size(); ++iNewPart) {
                except(newDomains[iNewPart], domainsWithId[iPrevious].domain, exceptedDomains);
            }
            newDomains.swap(exceptedDomains);
        }
        for (pluint iNew = 0; iNew < newDomains.size(); ++iNew) {
            nonOverlapingBlocks.push_back(
                DomainAndId3D(newDomains[iNew], domainsWithId[iDomain].id));
        }
    }
    return nonOverlapingBlocks;
}

void intersectDomainsAndIds(
    std::vector<std::vector<DomainAndId3D> > const &domainsWithId, std::vector<Box3D> &finalDomains,
    std::vector<std::vector<plint> > &finalIds)
{
    if (domainsWithId.empty()) {
        return;
    }

    // Copy the first collection without modification, as a starting
    //   point for future intersections.
    for (pluint iDomain = 0; iDomain < domainsWithId[0].size(); ++iDomain) {
        finalDomains.push_back(domainsWithId[0][iDomain].domain);
        finalIds.resize(finalDomains.size());
        finalIds.back().push_back(domainsWithId[0][iDomain].id);
    }

    // Then, intersect with all following collections.
    for (pluint iCollection = 1; iCollection < domainsWithId.size(); ++iCollection) {
        std::vector<Box3D> nextDomains;
        std::vector<std::vector<plint> > nextIds;
        // For each domain of the next collection ...
        for (pluint iNewDomain = 0; iNewDomain < domainsWithId[iCollection].size(); ++iNewDomain) {
            // ... and for each domain found so far ...
            for (pluint iOldDomain = 0; iOldDomain < finalDomains.size(); ++iOldDomain) {
                Box3D intersection;
                // ... compute the common intersections ...
                if (intersect(
                        domainsWithId[iCollection][iNewDomain].domain, finalDomains[iOldDomain],
                        intersection))
                {
                    // ... and store them.
                    nextDomains.push_back(intersection);
                    // As for the IDs, we need not only the ID of the newest block, but those of all
                    //   the previous blocks as well.
                    nextIds.push_back(finalIds[iOldDomain]);
                    nextIds.back().push_back(domainsWithId[iCollection][iNewDomain].id);
                }
            }
        }
        // Replace the old collections of intersections and IDs by the newest version.
        finalDomains.swap(nextDomains);
        finalIds.swap(nextIds);
    }
}

}  // namespace plb
