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


#ifndef IBDATAWRITINGFUNCTIONALS3D_HH_LBDEM
#define IBDATAWRITINGFUNCTIONALS3D_HH_LBDEM

#include "utils.h"

namespace plb {
  /* 
     implementation of getScalarQuantityFromDynamicsFunctional
  */

  template<typename T1, template<typename U> class Descriptor, typename T2>
  void GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>::
  process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
          ScalarField3D<T2>& data) 
  {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          T1 val = 0;
        
          Cell<T1,Descriptor>& cell = lattice.get(iX,iY,iZ);

          IBdynamicsParticleData<T1,Descriptor>* particleData =
            getParticleDataFromCell<T1,Descriptor>(cell);

          if(particleData){
            switch(which){
            case SolidFraction:
              val = (T2) particleData->solidFraction;
              break;
            case ParticleId:
              val = (T2) particleData->partId;
              break;
          
            }
          }
          data.get(iX,iY,iZ) = val;
        }
      }
    }
  }
  template<typename T1, template<typename U> class Descriptor, typename T2>
  GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>*      
  GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>::clone() const
  {
    return new GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>(*this);
  }
  
  template<typename T1, template<typename U> class Descriptor, typename T2>
  void GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>::
  getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
  }

  /* 
     implementation of getVectorQuantityFromDynamicsFunctional
  */

  template<typename T1, template<typename U> class Descriptor, typename T2, int nDim>
  void GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>::
  process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
          TensorField3D<T2,nDim>& data) 
  {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
            
          Array<T2,nDim> val;
          val.resetToZero();

          Cell<T1,Descriptor>& cell = lattice.get(iX,iY,iZ);

          IBdynamicsParticleData<T1,Descriptor>* particleData =
            getParticleDataFromCell<T1,Descriptor>(cell);

          if(particleData){
            switch(which){
            case ParticleVelocity:
              val = particleData->uPart;
              break;
            case HydrodynamicForce:
              val = particleData->hydrodynamicForce;
              break;
            
            }
          }
          data.get(iX,iY,iZ) = val;
        }
      }
    }
  }

  template<typename T1, template<typename U> class Descriptor, typename T2, int nDim>
  GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>* 
  GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>::clone() const
  {
    return new GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>(*this);
  }

  template<typename T1, template<typename U> class Descriptor, typename T2, int nDim>
  void GetVectorQuantityFromDynamicsFunctional<T1,Descriptor,T2,nDim>::
  getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
  }


}; /* namespace plb */

#endif /* IBDATAWRITINGFUNCTIONALS3D_HH_LBDEM */
