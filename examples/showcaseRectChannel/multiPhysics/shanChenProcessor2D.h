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

/* The original version of this file was written by Orestis Malaspinas
 * and Andrea Parmigiani.
 */

#ifndef SHAN_CHEN_PROCESSOR_2D_H
#define SHAN_CHEN_PROCESSOR_2D_H

#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "core/globalDefs.h"
#include "multiPhysics/interparticlePotential.h"

namespace plb {

/// Shan-Chen coupling for multi-component flow with or without external force
template <typename T, template <typename U> class Descriptor>
class ShanChenMultiComponentProcessor2D : public LatticeBoxProcessingFunctional2D<T, Descriptor> {
public:
    /// With these constructors, space- and time-dependent values of the
    ///   relaxation parameters omega are accounted for.
    ShanChenMultiComponentProcessor2D(T G_);
    ShanChenMultiComponentProcessor2D(std::vector<std::vector<T> > const &speciesG_);
    /// With these constructors, the values of the relaxation parameters omega are
    ///   taken to be species-dependent, but not space- or time-dependent. Their
    ///   value is imposed in the constructor.
    ShanChenMultiComponentProcessor2D(T G_, std::vector<T> const &imposedOmega_);
    ShanChenMultiComponentProcessor2D(
        std::vector<std::vector<T> > const &speciesG_, std::vector<T> const &imposedOmega_);
    virtual void process(Box2D domain, std::vector<BlockLattice2D<T, Descriptor> *> lattices);
    virtual ShanChenMultiComponentProcessor2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T G;
    std::vector<T> speciesG;
    std::vector<T> imposedOmega;
};

/// Shan-Chen coupling for single-component flow with or without external force
template <typename T, template <typename U> class Descriptor>
class ShanChenSingleComponentProcessor2D : public BoxProcessingFunctional2D_L<T, Descriptor> {
public:
    ShanChenSingleComponentProcessor2D(T G_, interparticlePotential::PsiFunction<T> *Psi_);
    virtual ~ShanChenSingleComponentProcessor2D();
    ShanChenSingleComponentProcessor2D(
        ShanChenSingleComponentProcessor2D<T, Descriptor> const &rhs);
    ShanChenSingleComponentProcessor2D &operator=(
        ShanChenSingleComponentProcessor2D<T, Descriptor> const &rhs);
    virtual void process(Box2D domain, BlockLattice2D<T, Descriptor> &lattice);
    virtual ShanChenSingleComponentProcessor2D<T, Descriptor> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;

private:
    T G;
    interparticlePotential::PsiFunction<T> *Psi;
};

}  // namespace plb

#endif  // SHAN_CHEN_LATTICES_2D_H
