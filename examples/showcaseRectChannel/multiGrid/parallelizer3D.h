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

#ifndef PARALLELIZER_3D_H
#define PARALLELIZER_3D_H

#include "core/geometry3D.h"

namespace plb {

class Parallelizer3D {
public:
    virtual ~Parallelizer3D() { }

    /// Method that will choose the parallelizing approach
    virtual void parallelize() = 0;

    virtual std::vector<std::vector<Box3D> > &getRecomputedBlocks()
    {
        return recomputedBlocks;
    }
    virtual std::vector<std::vector<Box3D> > const &getRecomputedBlocks() const
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

    /// Computes the "cost" of the domain delimited by b
    virtual plint computeCost(std::vector<std::vector<Box3D> > const &originalBlocks, Box3D b);

    virtual Parallelizer3D *clone() = 0;

    virtual void parallelizeLevel(
        plint whichLevel, std::vector<std::vector<Box3D> > const &originalBlocks,
        std::vector<Box3D> const &parallelRegions, std::vector<plint> const &regionIDs);

protected:
    // the results produced by the parallelizer
    std::vector<std::vector<Box3D> > recomputedBlocks;
    std::vector<std::vector<plint> > finalMpiDistribution;
};

/// Parallelize by cubes
class ParallellizeByCubes3D : public Parallelizer3D {
public:
    ParallellizeByCubes3D(
        std::vector<std::vector<Box3D> > const &originalBlocks_, Box3D finestBoundingBox_,
        plint xTiles_, plint yTiles_, plint zTiles_, bool optimiseDivision = false);

    virtual ~ParallellizeByCubes3D() { }

    /// Compute the new distribution of the blocks in the management
    virtual void parallelize();

    virtual Parallelizer3D *clone()
    {
        return new ParallellizeByCubes3D(originalBlocks, finestBoundingBox, xTiles, yTiles, zTiles);
    }

private:
    /// Create a regular division of the finest Box3D in nx x ny squares
    void computeFinestDivision();
    void computeDivision();

private:
    std::vector<std::vector<Box3D> > const &originalBlocks;
    Box3D finestBoundingBox;

    plint processorNumber;
    plint xTiles;
    plint yTiles;
    plint zTiles;

    // division the finest level in a number or regular squares
    std::vector<Box3D> finestDivision;
    // distribution w.r.t. finest division (same structure)
    std::vector<plint> mpiDistribution;
};

/// Parallelize by making an approximative load balancing
// class ParallellizeByLoadBalancing3D : public Parallelizer3D {
//     public:
//         ParallellizeByLoadBalancing3D(std::vector<std::vector<Box3D> > const& originalBlocks_,
//                                       Box3D finestBoundingBox_);

//        virtual ~ParallellizeByLoadBalancing3D(){}

//        /// Compute the new distribution of the blocks in the management
//        virtual void parallelize();

//        virtual Parallelizer3D* clone(){
//            return new ParallellizeByLoadBalancing3D(
//                        originalBlocks,finestBoundingBox);
//        }

// private:
//         void splitBlockInTwo(Box3D originalBox, Box3D& box1, Box3D& box2);

//    private:
//        std::vector<std::vector<Box3D> > const& originalBlocks;
//        Box3D finestBoundingBox;

//        plint processorNumber;
//};

}  // namespace plb

#endif  // PARALLELIZER_3D_H
