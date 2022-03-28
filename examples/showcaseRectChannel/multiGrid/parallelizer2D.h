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

#ifndef PARALLELIZER_2D_H
#define PARALLELIZER_2D_H

#include "core/geometry2D.h"

namespace plb {

class Parallelizer2D {
public:
    virtual ~Parallelizer2D() { }

    /// Method that will choose the parallelizing approach
    virtual void parallelize() = 0;

    virtual std::vector<std::vector<Box2D> > &getRecomputedBlocks()
    {
        return recomputedBlocks;
    }
    virtual std::vector<std::vector<Box2D> > const &getRecomputedBlocks() const
    {
        return recomputedBlocks;
    }

    virtual std::vector<std::vector<plint> > &getMpiDistribution()
    {
        return finalMpiDistribution;
    }
    virtual std::vector<std::vector<plint> > const &getMpiDistribution() const
    {
        return finalMpiDistribution;
    }

    virtual Parallelizer2D *clone() = 0;

    virtual void parallelizeLevel(
        plint whichLevel, std::vector<std::vector<Box2D> > const &originalBlocks,
        std::vector<Box2D> const &parallelRegions, std::vector<plint> const &regionIDs);

    /// Computes the "cost" of the domain delimited by b
    virtual plint computeCost(std::vector<std::vector<Box2D> > const &originalBlocks, Box2D b);

protected:
    // the results produced by the parallelizer
    std::vector<std::vector<Box2D> > recomputedBlocks;
    std::vector<std::vector<plint> > finalMpiDistribution;
};

/// Parallelize by making an approximative load balancing
class ParallellizeBySquares2D : public Parallelizer2D {
public:
    virtual ~ParallellizeBySquares2D() { }

    ParallellizeBySquares2D(
        std::vector<std::vector<Box2D> > const &originalBlocks_, Box2D finestBoundingBox_,
        plint xTiles, plint yTiles);
    /// Compute the new distribution of the blocks in the management
    virtual void parallelize();

    virtual Parallelizer2D *clone()
    {
        return new ParallellizeBySquares2D(originalBlocks, finestBoundingBox, xTiles, yTiles);
    }

private:
    /// Create a regular division of the finest Box2D in nx x ny squares
    void computeFinestDivision(plint xTiles, plint yTiles);

private:
    std::vector<std::vector<Box2D> > const &originalBlocks;
    Box2D finestBoundingBox;

    plint processorNumber;
    plint xTiles;
    plint yTiles;

    // division the finest level in a number or regular squares
    std::vector<Box2D> finestDivision;
    // distribution w.r.t. finest division (same structure)
    std::vector<plint> mpiDistribution;
};

}  // namespace plb

#endif  // PARALLELIZER_2D_H
