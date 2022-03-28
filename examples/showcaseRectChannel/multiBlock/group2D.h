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
 * 2D groups -- header file.
 */
#ifndef GROUP_2D_H
#define GROUP_2D_H

#include <map>
#include <string>
#include <vector>

#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiContainerBlock2D.h"
#include "multiBlock/multiDataField2D.h"
#include "particles/multiParticleField2D.h"

namespace plb {

class Group2D {
public:
    /// Attention: use carefully! There's a default constructor, for efficiency
    /// reasons. But the first thing you must do, before anything else, is to
    /// add a multi-block to the group.
    Group2D() { }
    Group2D(MultiBlock2D *block, std::string name = "");
    ~Group2D();
    Group2D(Group2D const &rhs);
    Group2D &operator=(Group2D const &rhs);
    void swap(Group2D &rhs);
    plint getNumBlocks() const;
    std::string getName(plint id) const;
    plint add(MultiBlock2D *block, std::string name = "");
    MultiBlock2D &get(plint id);
    MultiBlock2D &get(std::string name);
    bool hasBlock(std::string name) const;
    Box2D getBoundingBox() const;
    MultiBlockManagement2D const &getMultiBlockManagement() const;
    /// Replace the multi-block management of all multi-blocks in this group.
    void replaceManagement(MultiBlockManagement2D const &management);
    /// Replace the block with id "id", and re-assign the multi-block management
    /// of all other blocks accordingly.
    void replace(plint id, MultiBlock2D *block);
    /// Replace the block with name "name", and re-assign the multi-block management
    /// of all other blocks accordingly.
    void replace(std::string name, MultiBlock2D *block);

    template <typename T>
    plint generateScalar(std::string name = "", plint envelopeWidth = 1, plint gridLevel = 0);
    template <typename T>
    plint generateNTensor(
        std::string name = "", plint nDim = 1, plint envelopeWidth = 1, plint gridLevel = 0);
    template <typename T, int nDim>
    plint generateTensor(std::string name = "", plint envelopeWidth = 1, plint gridLevel = 0);
    template <typename T, template <typename U> class Descriptor>
    plint generateLattice(std::string name = "", plint envelopeWidth = 1, plint gridLevel = 0);
    template <class ParticleFieldT>
    plint generateParticleField(
        std::string name = "", plint envelopeWidth = 1, plint gridLevel = 0);
    plint generateContainer(std::string name = "", plint envelopeWidth = 1, plint gridLevel = 0);

    template <typename T>
    MultiScalarField2D<T> &getScalar(plint id);
    template <typename T>
    MultiScalarField2D<T> &getScalar(std::string name);
    template <typename T>
    MultiNTensorField2D<T> &getNTensor(plint id);
    template <typename T>
    MultiNTensorField2D<T> &getNTensor(std::string name);
    template <typename T, int nDim>
    MultiTensorField2D<T, nDim> &getTensor(plint id);
    template <typename T, int nDim>
    MultiTensorField2D<T, nDim> &getTensor(std::string name);
    template <typename T, template <typename U> class Descriptor>
    MultiBlockLattice2D<T, Descriptor> &getLattice(plint id);
    template <typename T, template <typename U> class Descriptor>
    MultiBlockLattice2D<T, Descriptor> &getLattice(std::string name);
    template <class ParticleFieldT>
    MultiParticleField2D<ParticleFieldT> &getParticleField(plint id);
    template <class ParticleFieldT>
    MultiParticleField2D<ParticleFieldT> &getParticleField(std::string name);
    MultiContainerBlock2D &getContainer(plint id);
    MultiContainerBlock2D &getContainer(std::string name);

private:
    std::vector<MultiBlock2D *> blocks;
    std::map<std::string, plint> ids;
    std::map<plint, std::string> names;
};

}  // namespace plb

#endif  // GROUP_2D_H
