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
 * 3D groups -- header file.
 */
#ifndef GROUP_3D_H
#define GROUP_3D_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiDataField3D.h"
#include "particles/multiParticleField3D.h"

namespace plb {

class Group3D {
public:
    /// Attention: use carefully! There's a default constructor, for efficiency
    /// reasons. But the first thing you must do, before anything else, is to
    /// add a multi-block to the group.
    Group3D() { }
    Group3D(MultiBlock3D *block, std::string name = "");
    ~Group3D();
    Group3D(Group3D const &rhs);
    Group3D &operator=(Group3D const &rhs);
    void swap(Group3D &rhs);
    std::string getName(plint id) const;
    plint add(MultiBlock3D *block, std::string name = "");
    plint getNumBlocks() const;
    MultiBlock3D &get(plint id);
    MultiBlock3D &get(std::string name);
    bool hasBlock(std::string name) const;
    Box3D getBoundingBox() const;
    MultiBlockManagement3D const &getMultiBlockManagement() const;
    /// Replace the multi-block management of all multi-blocks in this group.
    void replaceManagement(MultiBlockManagement3D const &management);
    /// Replace the block with id "id", and re-assign the multi-block management
    /// of all other blocks accordingly.
    void replace(plint id, MultiBlock3D *block);
    /// Replace the block with name "name", and re-assign the multi-block management
    /// of all other blocks accordingly.
    void replace(std::string name, MultiBlock3D *block);

    template <typename T>
    plint generateScalar(std::string name = "", plint envelopeWidth = 1, plint gridLevel = 0);
    template <typename T>
    plint generateNTensor(
        std::string name = "", plint nDim = 1, plint envelopeWidth = 1, plint gridLevel = 0);
    template <typename T>
    plint generateNTensor(
        std::string name, plint nDim, std::vector<Box3D> const &bulks, plint envelopeWidth,
        plint gridLevel, bool cropBoundingBox);
    template <typename T, int nDim>
    plint generateTensor(std::string name = "", plint envelopeWidth = 1, plint gridLevel = 0);
    template <typename T, template <typename U> class Descriptor>
    plint generateLattice(std::string name = "", plint envelopeWidth = 1, plint gridLevel = 0);
    template <typename T, template <typename U> class Descriptor>
    plint generateDenseParticles(
        std::string name = "", plint envelopeWidth = 1, plint gridLevel = 0);
    plint generateContainer(std::string name = "", plint envelopeWidth = 1, plint gridLevel = 0);
    plint generateContainerWithData(
        std::string name = "", ContainerBlockData *data = 0, plint envelopeWidth = 1,
        plint gridLevel = 0);

    template <typename T>
    MultiScalarField3D<T> &getScalar(plint id);
    template <typename T>
    MultiScalarField3D<T> &getScalar(std::string name);
    template <typename T>
    MultiNTensorField3D<T> &getNTensor(plint id);
    template <typename T>
    MultiNTensorField3D<T> &getNTensor(std::string name);
    template <typename T, int nDim>
    MultiTensorField3D<T, nDim> &getTensor(plint id);
    template <typename T, int nDim>
    MultiTensorField3D<T, nDim> &getTensor(std::string name);
    template <typename T, template <typename U> class Descriptor>
    MultiBlockLattice3D<T, Descriptor> &getLattice(plint id);
    template <typename T, template <typename U> class Descriptor>
    MultiBlockLattice3D<T, Descriptor> &getLattice(std::string name);
    template <typename T, template <typename U> class Descriptor>
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &getDenseParticles(plint id);
    template <typename T, template <typename U> class Descriptor>
    MultiParticleField3D<DenseParticleField3D<T, Descriptor> > &getDenseParticles(std::string name);
    MultiContainerBlock3D &getContainer(plint id);
    MultiContainerBlock3D &getContainer(std::string name);

private:
    plint addNoCheck(MultiBlock3D *block, std::string name = "");

private:
    std::vector<MultiBlock3D *> blocks;
    std::map<std::string, plint> ids;
    std::map<plint, std::string> names;
};

template <typename T, typename TConv>
void addTransform(
    Group3D &group, std::unique_ptr<MultiScalarField3D<T> > field, std::string name,
    TConv scalingFactor, TConv additiveOffset = (TConv)0);

template <typename T>
void addTransform(
    Group3D &group, std::unique_ptr<MultiScalarField3D<T> > field, std::string name,
    T scalingFactor, T additiveOffset = (T)0);

template <typename T, typename TConv, int nDim>
void addTransform(
    Group3D &group, std::unique_ptr<MultiTensorField3D<T, nDim> > field, std::string name,
    TConv scalingFactor);

template <typename T, int nDim>
void addTransform(
    Group3D &group, std::unique_ptr<MultiTensorField3D<T, nDim> > field, std::string name,
    T scalingFactor);

template <typename T, typename TConv>
void addTransform(
    Group3D &group, MultiScalarField3D<T> &field, std::string name, TConv scalingFactor,
    TConv additiveOffset = (TConv)0);

template <typename T>
void addTransform(
    Group3D &group, MultiScalarField3D<T> &field, std::string name, T scalingFactor,
    T additiveOffset = (T)0);

template <typename T, typename TConv, int nDim>
void addTransform(
    Group3D &group, MultiTensorField3D<T, nDim> &field, std::string name, TConv scalingFactor);

template <typename T, int nDim>
void addTransform(
    Group3D &group, MultiTensorField3D<T, nDim> &field, std::string name, T scalingFactor);

}  // namespace plb

#endif  // GROUP_3D_H
