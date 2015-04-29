/*
 * This file is part of the LBDEMcoupling software.
 *
 * LBDEMcoupling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2014 Johannes Kepler University Linz
 *
 * Author: Philippe Seil (philippe.seil@jku.at)
 */

/*
 * functionals used to obtain particle data that can be used for post processing
 */

#ifndef IBDATAWRITINGFUNCTIONALS3D_H_LBDEM
#define IBDATAWRITINGFUNCTIONALS3D_H_LBDEM

namespace plb {

  enum IBscalarQuantity { SolidFraction, ParticleId };
  enum IBvectorQuantity { ParticleVelocity, HydrodynamicForce };

  template<typename T1, template<typename U> class Descriptor, typename T2>
  struct GetScalarQuantityFromDynamicsFunctional : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2> {
  public:

    IBscalarQuantity which;
    GetScalarQuantityFromDynamicsFunctional(IBscalarQuantity const which_)
      : which(which_) {}
    
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                         ScalarField3D<T2>& data);

    GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>* clone() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    
  };

  template<typename T1, template<typename U> class Descriptor, typename T2, int nDim>
  struct GetVectorQuantityFromDynamicsFunctional : public BoxProcessingFunctional3D_LT<T1,Descriptor,T2,nDim> {
  public:

    IBvectorQuantity which;
    GetVectorQuantityFromDynamicsFunctional(IBvectorQuantity const which_)
      : which(which_) {}
    
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                         TensorField3D<T2,nDim>& data);

    GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>* clone() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
  };

}; /* namespace plb */

#include "ibDataWritingFunctionals3D.hh"

#endif /* IBDATAWRITINGFUNCTIONALS3D_H_LBDEM */
